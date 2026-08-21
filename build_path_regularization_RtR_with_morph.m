function [RtR, diagnostics] = build_path_regularization_RtR_with_morph(imdl, thr_factor, reduce_factor, morph_thresh, varargin)
% BUILD_PATH_REGULARIZATION_RTR_WITH_MORPH  Path detection + morphological reconstruction
% RtR = build_path_regularization_RtR_with_morph(imdl, thr_factor, reduce_factor, morph_thresh)
%
% See Copy_3_of_cc33.m for usage. This function is factored out so it is
% visible on the MATLAB path when called inside inv_solve.
if nargin < 2 || isempty(thr_factor), thr_factor = 0.5; end
if nargin < 3 || isempty(reduce_factor), reduce_factor = 0.2; end
if nargin < 4 || isempty(morph_thresh), morph_thresh = 0.3; end

% optional params
% numerical / behaviour options
% - ridge: small diagonal term added to RtR to avoid near-singularity
% - blend_uniform: blend fraction toward the mean weight to avoid isolated near-zero rows
opts = struct('p',2,'inner_iters',3,'alpha',0.7,'avg_iters',1,'lambda',1,'do_opening',false,'do_closing',true, ...
              'w_floor',0.2,'w_ceil',1.5,'ridge',1e-4,'blend_uniform',0.08, ...
              'use_percentile',true,'hthr_pct',85,'lthr_pct',60,'max_path_frac',0.15,'max_path_abs',300, ...
              'path_mode','reduce', ... % 'reduce' 或 'boost'
              'boost_factor',0.5, ...  % path_mode='boost' 时使用
              'combine_mode','plain', ... % 'plain' RtR = (W L)^T (W L); 'blend' 将与 TV 混合
              'blend_alpha',0.5, ... % combine_mode='blend' 的混合系数
              'normalize_mode','trace', ... % 对 Lw 进行尺度归一：'none' | 'trace' | 'fro'
              'use_midline_band', true, ... % 优先在中轴带内选端点以得到中央连通路径
              'midline_band', 0.25); % |y| < band 的元素优先用于端点选择
if ~isempty(varargin)
    user = varargin{1}; fn = fieldnames(user);
    for k=1:numel(fn), opts.(fn{k}) = user.(fn{k}); end
end

% 获取元素初始值
if isfield(imdl, 'elem_data')
    base_elem_data = imdl.elem_data;
elseif isfield(imdl, 'jacobian_bkgnd') && isfield(imdl.jacobian_bkgnd, 'value')
    base_elem_data = imdl.jacobian_bkgnd.value * ones(size(prior_TV(imdl),2),1);
else
    L = prior_TV(imdl); RtR = L' * L; return;
end

% 基础算子
L = prior_TV(imdl);
EInc = spones(abs(L));
G = size(L,1);
Ne = size(L, 2);  % 元素数量
edges_per_elem = full(max(1, sum(EInc,1)'));

% 元素邻接矩阵
A = spones(EInc' * EInc);
A = A - spdiags(diag(A),0,size(A,1),size(A,2));
deg = full(sum(A,2));

% ========== 按论文：使用梯度算子检测边缘 ==========
% 论文: "计算梯度幅值场 g = |L * σ|"
% g_edge: 每条边上的梯度（相邻元素电导率之差的绝对值）
g_edge = abs(L * base_elem_data);

% 将边梯度聚合到元素域
g_elem = (EInc' * g_edge) ./ edges_per_elem;

% 邻接平滑
for t = 1:opts.avg_iters
    g_elem = (A * g_elem + opts.lambda * g_elem) ./ max(1, (deg + opts.lambda));
end

% 形态学开闭运算（论文: "执行形态学开闭运算以滤除噪声"）
if opts.do_opening
    g_elem = element_graph_erosion(g_elem, A);
    g_elem = element_graph_dilation(g_elem, A);
end
if opts.do_closing
    g_elem = element_graph_dilation(g_elem, A);
    g_elem = element_graph_erosion(g_elem, A);
end

% 双阈值分割（论文: "自适应计算高低阈值"）
if isfield(opts,'use_percentile') && opts.use_percentile
    hthr = prctile(g_elem, opts.hthr_pct);
    lthr = prctile(g_elem, opts.lthr_pct);
else
    hthr = mean(g_elem) + thr_factor * std(g_elem);
    lthr = min(g_elem) + morph_thresh * (max(g_elem)-min(g_elem));
end

% 生成 marker 和 mask（论文: "生成标记图像 M 与掩模图像 F"）
marker = g_elem > hthr;  % 高梯度核心区域
mask = g_elem > lthr;    % 较低阈值的候选区域

if sum(marker) < 1
    RtR = L' * L; 
    if nargout > 1
        diagnostics = struct('path_nodes', [], 'path_length', 0, 'g_elem', g_elem, ...
            'g_elem2', g_elem, 'mask2', false(Ne,1), 'recon', false(Ne,1), ...
            'hthr', hthr, 'max_path_len', 0, 'w', ones(G,1));
    end
    return;
end

% 测地线形态学重建（论文: "执行测地线形态学重建 R = ρ_F(M)"）
recon = graph_morphological_reconstruction(marker, mask, A);

% g_elem2 用于诊断
g_elem2 = g_elem;

% mask2 = 形态学重建结果（高梯度连通区域 = 边缘/骨架）
mask2 = recon;

if sum(mask2) < 2
    RtR = L' * L; 
    if nargout > 1
        diagnostics = struct('path_nodes', [], 'path_length', 0, 'g_elem', g_elem, ...
            'g_elem2', g_elem2, 'mask2', mask2, 'recon', recon, ...
            'hthr', hthr, 'max_path_len', 0, 'w', ones(G,1));
    end
    return;
end

% ========== 接地网区域选择 ==========
% 对于接地网成像，应该使用**所有满足阈值的区域**，而不仅仅是最大连通分量
% 因为接地网的边界在 pilot 重建中可能是断开的多个区域
% 选项：use_all_components = true 时使用所有区域，false 时只用最大连通分量
use_all_components = true;

idx = find(mask2);

if use_all_components
    % 使用所有满足阈值的元素（接地网的多个边界区域）
    path_nodes = idx;
else
    % 原逻辑：只提取最大连通分量
    A_sub = A(idx, idx);
    Gsub = graph(A_sub);
    if numnodes(Gsub) == 0
        RtR = L' * L; return;
    end
    
    cc = conncomp(Gsub);
    unique_cc = unique(cc);
    counts = arrayfun(@(c) sum(cc==c), unique_cc);
    [~, imax] = max(counts);
    nodes_sub = find(cc == unique_cc(imax));
    
    if numel(nodes_sub) < 2
        RtR = L' * L; return;
    end
    
    path_nodes = idx(nodes_sub);
end

% 可选：区域大小限制
max_len = min(opts.max_path_abs, max(10, ceil(opts.max_path_frac * Ne)));
if numel(path_nodes) > max_len
    % 按梯度幅值排序，保留最高的 max_len 个元素
    [~, sort_idx] = sort(g_elem(path_nodes), 'descend');
    path_nodes = path_nodes(sort_idx(1:max_len));
end

% 将路径映射回边
onpath = false(G,1);
for rr = 1:G
    elems = find(L(rr,:) ~= 0);
    if any(ismember(elems, path_nodes))
        onpath(rr) = true;
    end
end

% IRLS 风格的内循环：基于加权算子迭代更新边权
w = ones(G,1);
% set a tau for weight function based on edge gradient
g_edge_for_tau = abs(L * base_elem_data);
tau = max(median(g_edge_for_tau), eps);
for iter = 1:opts.inner_iters
    Wtmp = spdiags(w,0,G,G);
    Lw_tmp = Wtmp * L;
    g_edge_tmp = abs(Lw_tmp * base_elem_data);
    % 元素域聚合
    g_elem_tmp = (EInc' * g_edge_tmp) ./ edges_per_elem;
    for t = 1:opts.avg_iters
        g_elem_tmp = (A * g_elem_tmp + opts.lambda * g_elem_tmp) ./ max(1, (deg + opts.lambda));
    end
    if opts.do_opening
        g_elem_tmp = element_graph_erosion(g_elem_tmp, A);
        g_elem_tmp = element_graph_dilation(g_elem_tmp, A);
    end
    if opts.do_closing
        g_elem_tmp = element_graph_dilation(g_elem_tmp, A);
        g_elem_tmp = element_graph_erosion(g_elem_tmp, A);
    end
    % 回映到边域
    g_edge_ref = zeros(G,1);
    for rr = 1:G
        elems = find(L(rr,:) ~= 0);
        if isempty(elems)
            g_edge_ref(rr) = 0;
        else
            g_edge_ref(rr) = mean(g_elem_tmp(elems));
        end
    end
    w_new = 1 ./ (1 + (g_edge_ref ./ max(tau,eps)).^opts.p);
    w = opts.alpha * w + (1-opts.alpha) * w_new;
end

% 基于图距离生成连续路径衰减权重（比二值 onpath 更平滑）
% 1) 在元素图上计算到 path_nodes 的最短路径距离
Gfull = graph(A);
% distances accepts vector of source nodes; returns matrix (#sources x #nodes)
Dmat = distances(Gfull, path_nodes);
if isempty(Dmat)
    dist_elem = inf(size(A,1),1);
else
    dist_elem = min(Dmat, [], 1)'; % Ne x 1
end

% 2) 将元素距离映射回边域（取两端元素的平均距离）
dist_edge = zeros(G,1);
for rr = 1:G
    elems = find(L(rr,:) ~= 0);
    if isempty(elems)
        dist_edge(rr) = inf;
    else
        dist_edge(rr) = mean(dist_elem(elems));
    end
end

% 3) 路径相关的权重缩放
% sigma 控制衰减/提升宽度（以元素单位）
sigma = max(1, round(max(1, mean(dist_elem(dist_elem<inf))/2))); % 自适应默认
if ~isfinite(sigma) || sigma<=0, sigma = 3; end

% ========== 修正：权重逻辑 ==========
% 目标：接地网区域内 TV 惩罚小（允许高电导率/锐利边缘）
%       接地网区域外 TV 惩罚大（保持背景平滑）
%
% 权重 w 小 → Lw = W*L 小 → RtR = Lw'*Lw 小 → TV 惩罚小
% 所以：接地网区域内应该权重小，区域外权重大
%
switch lower(string(opts.path_mode))
    case "boost"
        % 提升路径外侧惩罚（接地网内权重=1，外侧权重增加）
        bf = opts.boost_factor; if ~isfinite(bf) || bf<0, bf = 0.5; end
        path_scale = 1 + bf * (1 - exp( - (dist_edge./sigma).^2 ));
    otherwise
        % 默认 'reduce'：降低接地网区域内的权重
        %   在接地网内 (dist=0): scale = reduce_factor（小）→ TV惩罚小
        %   远离接地网 (dist大): scale → 1（大）→ TV惩罚大
        path_scale = 1 - (1 - reduce_factor) * exp( - (dist_edge./sigma).^2 );
end

% 4) combine IRLS weight w with path_scale
% keep a copy of the IRLS-only weight for diagnostics
w_irls = w;
w = w .* path_scale;

% 5) 邻接平滑权重（在边域上，缓和孤立突变）
% 构造边邻接：若两条边共享单元则相邻（G x G）
A_edge = spones(L * L');
deg_edge = full(sum(A_edge,2));
for t = 1:2
    w = (A_edge * w + w) ./ max(1, (deg_edge + 1));
end

% blend slightly toward the mean weight to avoid isolated tiny rows/eigenvalues
if isfinite(opts.blend_uniform) && opts.blend_uniform > 0
    wm = mean(w(isfinite(w))); if isempty(wm) || ~isfinite(wm), wm = 1; end
    w = opts.blend_uniform * wm + (1-opts.blend_uniform) * w;
end

% guard against non-finite values
w(~isfinite(w)) = min(1, max(0, nanmedian(w)));

% 截断并返回 RtR
w = (w - min(w)) / max(eps, (max(w)-min(w))); % [0,1]
w = opts.w_floor + w .* (opts.w_ceil - opts.w_floor);
W = spdiags(w,0,G,G);
Lw = W * L;
% 可选尺度归一：让加权算子的整体尺度接近原始 L，从而 blend_alpha 的解释更稳定
scale_info = struct('mode', opts.normalize_mode, 's', 1);
switch lower(string(opts.normalize_mode))
    case 'trace'
        tL = full(trace(L' * L)); tW = full(trace((Lw' * Lw)));
        if isfinite(tL) && isfinite(tW) && tW>0
            s = sqrt(tL / tW); Lw = s * Lw; scale_info.s = s; end
    case 'fro'
        fL = norm(L,'fro'); fW = norm(Lw,'fro');
        if isfinite(fL) && isfinite(fW) && fW>0
            s = fL / fW; Lw = s * Lw; scale_info.s = s; end
    otherwise
        scale_info.s = 1;
end
switch lower(string(opts.combine_mode))
    case "blend"
        % 与纯 TV 先验混合，降低偏差： RtR = a*(L^T L) + (1-a)*(Lw^T Lw)
        a = opts.blend_alpha; if ~isfinite(a) || a<0, a=0.5; elseif a>1, a=1; end
        RtR = a*(L' * L) + (1-a)*(Lw' * Lw) + opts.ridge * speye(size(L,2));
    otherwise
        RtR = Lw' * Lw + opts.ridge * speye(size(L,2));
end

% optional diagnostic output
if nargout > 1
    diagnostics = struct();
    diagnostics.w = w;
    diagnostics.w_irls = w_irls;      % IRLS-only weight (before path scaling)
    diagnostics.path_nodes = path_nodes;
    diagnostics.mask2 = mask2;
    diagnostics.recon = recon;
    diagnostics.g_elem = g_elem;
    diagnostics.g_elem2 = g_elem2;
    diagnostics.dist_edge = dist_edge;
    diagnostics.path_scale = path_scale;
    diagnostics.sigma = sigma;
    diagnostics.combine_mode = opts.combine_mode;
    diagnostics.blend_alpha = opts.blend_alpha;
    diagnostics.normalize_mode = opts.normalize_mode;
    diagnostics.scale_s = scale_info.s;
    % threshold diagnostics
    diagnostics.hthr = hthr;
    if isfield(opts,'use_percentile') && opts.use_percentile
        diagnostics.hthr_pct = opts.hthr_pct;
        diagnostics.lthr_pct = opts.lthr_pct;
    end
    diagnostics.max_path_len = max_len;
    diagnostics.path_length = numel(path_nodes);
else
    diagnostics = [];
end
end
