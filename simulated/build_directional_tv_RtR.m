function RtR = build_directional_tv_RtR(imdl, tau, beta, w_min, w_max)
% BUILD_DIRECTIONAL_TV_RTR  Directional-adaptive TV prior RtR builder (simplified)
% RtR = build_directional_tv_RtR(imdl, tau, beta, w_min, w_max)
%
% Parameters:
%  tau   - threshold for adaptive weighting (default 0.2)
%  beta  - directional sensitivity (default 0.8)
%  w_min - minimum weight (default 0.1)
%  w_max - maximum weight (default 1.0)

if nargin < 2 || isempty(tau), tau = 0.2; end
if nargin < 3 || isempty(beta), beta = 0.8; end
if nargin < 4 || isempty(w_min), w_min = 0.1; end
if nargin < 5 || isempty(w_max), w_max = 1.0; end

% 1) 获取基础元素数据
if isfield(imdl, 'elem_data')
    base_elem_data = imdl.elem_data;
elseif isfield(imdl, 'jacobian_bkgnd') && isfield(imdl.jacobian_bkgnd, 'value')
    base_elem_data = imdl.jacobian_bkgnd.value * ones(size(prior_TV(imdl),2),1);
else
    % 降级到标准TV
    L = prior_TV(imdl); 
    RtR = L' * L; 
    return;
end

% 2) 构建TV算子并计算梯度
L = prior_TV(imdl);                 % G x Ne (edges x elements)
G = size(L,1);  Ne = size(L,2);
g_edge = abs(L * base_elem_data);   % 每条边的梯度强度

% 3) 计算每个元素的平均梯度（简单聚合）
EInc = spones(abs(L));              % G x Ne 边-元素关联矩阵
edges_per_elem = full(max(1, sum(EInc,1)'));
g_elem = (EInc' * g_edge) ./ edges_per_elem; % Ne x 1

% 4) 提取元素中心坐标（用于方向估计）
centers = zeros(Ne, 2);
use_geom = false;
try
    nodes = imdl.fwd_model.nodes;
    elems = imdl.fwd_model.elems;
    if ~isempty(nodes) && ~isempty(elems) && size(nodes,2) >= 2
        for i = 1:Ne
            centers(i,:) = mean(nodes(elems(i,1:min(end,3)),1:2), 1);
        end
        use_geom = true;
    end
catch
    use_geom = false;
end

% 5) 估计每个元素的梯度主方向（简化版）
A = spones(EInc' * EInc);  % 元素邻接矩阵
A = A - spdiags(diag(A), 0, size(A,1), size(A,2));
dirs = zeros(Ne, 2);
for i = 1:Ne
    nb = find(A(i,:));
    if isempty(nb)
        dirs(i,:) = [0 0]; 
        continue;
    end
    gvec = [0 0];
    for kk = 1:numel(nb)
        j = nb(kk);
        gval = g_elem(j) - g_elem(i);
        if use_geom
            dv = centers(j,:) - centers(i,:);
        else
            dv = [j-i, 0];
        end
        if norm(dv) > eps
            gvec = gvec + gval * (dv / norm(dv));
        end
    end
    if norm(gvec) > eps
        dirs(i,:) = gvec / norm(gvec);
    end
end

% 6) 计算边权重（自适应 + 方向性）
w = ones(G, 1);
for rr = 1:G
    elems = find(L(rr,:) ~= 0);
    if isempty(elems)
        continue;
    end
    
    % 自适应权重（基于梯度强度）
    g_ref = mean(g_elem(elems));
    w_adaptive = 1 / (1 + (g_ref / max(tau, eps))^2);
    
    % 方向性权重
    if use_geom && numel(elems) == 2
        v = centers(elems(2),:) - centers(elems(1),:);
        ued = v / max(norm(v), eps);
    else
        ued = [1 0];
    end
    dir_avg = mean(dirs(elems,:), 1);
    if norm(dir_avg) > eps
        dir_avg = dir_avg / norm(dir_avg);
    end
    dotv = abs(dot(ued, dir_avg));
    if isnan(dotv), dotv = 0; end
    w_directional = 1 - beta * dotv;
    
    % 组合权重
    w(rr) = w_adaptive * w_directional;
end

% 诊断：归一化前的权重
w_raw_min = min(w); w_raw_max = max(w); w_raw_std = std(w);

% 7) 归一化并截断到 [w_min, w_max]
w_range = max(w) - min(w);
if w_range > eps
    w = (w - min(w)) / w_range;  % [0,1]
    w = w_min + w * (w_max - w_min);
else
    % 如果权重无变化，使用默认均匀权重
    w = ones(G,1) * mean([w_min, w_max]);
end

% 8) 构建加权 RtR 并打印诊断信息
W = spdiags(w, 0, G, G);
Lw = W * L;
RtR = Lw' * Lw;

% 诊断输出：权重统计（包括归一化前后）
fprintf('  Raw weights: min=%.4f, max=%.4f, std=%.4f | Final: min=%.4f, max=%.4f, mean=%.4f, std=%.4f\n', ...
    w_raw_min, w_raw_max, w_raw_std, min(w), max(w), mean(w), std(w));
end

