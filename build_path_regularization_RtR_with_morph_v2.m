function [RtR, diagnostics] = build_path_regularization_RtR_with_morph_v2(imdl, varargin)
% BUILD_PATH_REGULARIZATION_RTR_WITH_MORPH_V2  简化版形态学路径TV正则化
%
% 用法: RtR = build_path_regularization_RtR_with_morph_v2(imdl, opts)
%
% 参数 opts (struct):
%   ridge       - RtR对角项 (默认 1e-6)
%   w_floor     - 权重下限 (默认 0.25)
%   w_ceil      - 权重上限 (默认 1.0)
%   blend_alpha - 与标准TV混合 (默认 0.4)
%   hthr_pct    - 高阈值百分位 (默认 70)
%   lthr_pct    - 低阈值百分位 (默认 60)

%% 默认参数
opts = struct('ridge', 1e-10, 'w_floor', 0.25, 'w_ceil', 1.0, ...
              'blend_alpha', 0.4, 'hthr_pct', 50, 'lthr_pct', 30);
if ~isempty(varargin) && isstruct(varargin{end})
    user = varargin{end}; fn = fieldnames(user);
    for k = 1:numel(fn), opts.(fn{k}) = user.(fn{k}); end
end

%% 获取pilot电导率
if isfield(imdl, 'elem_data')
    sigma = imdl.elem_data;
else
    L = prior_TV(imdl); 
    RtR = L' * L; 
    diagnostics = struct('path_nodes',[],'path_length',0,'g_elem',[],'w',[]); 
    return;
end

%% 基础TV算子
L = prior_TV(imdl);
G = size(L, 1);   % 边数
Ne = size(L, 2);  % 元素数

% 元素邻接矩阵
EInc = spones(abs(L));
edges_per_elem = full(max(1, sum(EInc, 1)'));
A = spones(EInc' * EInc);
A = A - spdiags(diag(A), 0, Ne, Ne);

%% 1. 计算梯度场 g = |L*σ|
g_edge = abs(L * sigma);
g_elem = (EInc' * g_edge) ./ edges_per_elem;

% 形态学闭运算
g_elem = element_graph_dilation(g_elem, A);
g_elem = element_graph_erosion(g_elem, A);

%% 2. 双阈值 + 形态学重建
hthr = prctile(g_elem, opts.hthr_pct);
lthr = prctile(g_elem, opts.lthr_pct);
marker = g_elem > hthr;
mask = g_elem > lthr;

if sum(marker) < 1
    RtR = L' * L + opts.ridge * speye(Ne);
    diagnostics = struct('path_nodes',[],'path_length',0,'g_elem',g_elem,...
        'mask2',false(Ne,1),'hthr',hthr,'w',ones(G,1));
    return;
end

% 测地线形态学重建
recon = graph_morphological_reconstruction(marker, mask, A);
path_nodes = find(recon);  % 使用所有满足条件的区域

if numel(path_nodes) < 2
    RtR = L' * L + opts.ridge * speye(Ne);
    diagnostics = struct('path_nodes',[],'path_length',0,'g_elem',g_elem,...
        'mask2',recon,'hthr',hthr,'w',ones(G,1));
    return;
end

%% 3. 计算到骨架的距离
Gfull = graph(A);
Dmat = distances(Gfull, path_nodes);
dist_elem = min(Dmat, [], 1)';

% 映射到边域
dist_edge = zeros(G, 1);
for rr = 1:G
    elems = find(L(rr,:) ~= 0);
    if ~isempty(elems)
        dist_edge(rr) = mean(dist_elem(elems));
    else
        dist_edge(rr) = inf;
    end
end

%% 4. 计算权重
% 目标：骨架区域TV惩罚小 → 允许大梯度 → 边缘锐利
% 方法：骨架处权重小 → Lw行小 → ||Lw*σ||小 → 惩罚小
sigma_dist = max(1, mean(dist_elem(dist_elem < inf)) / 2);
if ~isfinite(sigma_dist), sigma_dist = 3; end

% w = w_floor + (w_ceil-w_floor) * (1 - exp(-d²/σ²))
% d=0（骨架上）: w = w_floor（小）→ TV惩罚小
% d大（远离骨架）: w → w_ceil（大）→ TV惩罚大
w = opts.w_floor + (opts.w_ceil - opts.w_floor) * (1 - exp(-(dist_edge / sigma_dist).^2));

% 边域平滑
A_edge = spones(L * L');
deg_edge = full(sum(A_edge, 2));
for t = 1:2
    w = (A_edge * w + w) ./ max(1, deg_edge + 1);
end

% 截断
w = max(opts.w_floor, min(opts.w_ceil, w));
w(~isfinite(w)) = opts.w_ceil;

%% 5. 构造RtR
% 方法：权重直接作用在惩罚矩阵对角上
W = spdiags(w, 0, G, G);
RtR_weighted = L' * W * L;

% --- Ridge (对角加载) ---
Ridge = opts.ridge * speye(Ne);
% ---------------------------------------------

% 混合（可选）
% blend_alpha 控制混合比例：
% 0.0 = 纯加权TV (RtR_weighted)
% 1.0 = 纯标准正则化 (Standard_Reg)
% 中间值 = 线性混合
a = opts.blend_alpha;

if a > 0
    % 定义标准正则化项 (Standard_Reg)
    % 选项 A: 标准 TV (L'L) - 保持边缘锐利，但可能阶梯化
    Standard_Reg = L' * L;
    
    % 选项 B: NOSER (加权 Laplacian) - 考虑灵敏度，中心惩罚大，边界小
    % Standard_Reg = prior_noser(imdl);
    
    % 选项 C: Tikhonov (I) - 平滑，抑制噪声，但模糊边缘
    % Standard_Reg = speye(Ne);
    
    % 选项 D: Laplacian (L_smooth' * L_smooth) - 平滑且连续
    % Standard_Reg = prior_laplace(imdl);
    
    % 直接混合 NOSER 和 weighted TV (量级相近，无需额外缩放)
    RtR = (1 - a) * RtR_weighted + a * Standard_Reg + Ridge;
else
    RtR = RtR_weighted + Ridge;
end

%% 诊断输出（合并字段，不覆盖先前诊断信息）
if nargout > 1
    if ~exist('diagnostics','var') || isempty(diagnostics)
        diagnostics = struct();
    end
    diagnostics.w = w;
    diagnostics.path_nodes = path_nodes;
    diagnostics.path_length = numel(path_nodes);
    diagnostics.mask2 = recon;
    diagnostics.g_elem = g_elem;
    diagnostics.hthr = hthr;
    diagnostics.lthr = lthr;
end
end
