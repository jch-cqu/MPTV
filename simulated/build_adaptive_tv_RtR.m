function RtR = build_adaptive_tv_RtR(imdl, tau, p, w_min, w_max)
% BUILD_ADAPTIVE_TV_RTR  Adaptive anisotropic TV RtR builder (standalone)
% RtR = build_adaptive_tv_RtR(imdl, tau, p, w_min, w_max)
%
% Parameters:
% tau  - scale for edge sensitivity
% p    - power controlling weight curve
% w_min, w_max - clipping range for weights

if nargin < 2 || isempty(tau), tau = 0.002; end
if nargin < 3 || isempty(p), p = 2.0; end
if nargin < 4 || isempty(w_min), w_min = 0.05; end
if nargin < 5 || isempty(w_max), w_max = 1.5; end

% obtain element baseline
if isfield(imdl, 'elem_data')
    base_elem_data = imdl.elem_data;
elseif isfield(imdl, 'jacobian_bkgnd') && isfield(imdl.jacobian_bkgnd, 'value')
    base_elem_data = imdl.jacobian_bkgnd.value * ones(size(prior_TV(imdl),2),1);
else
    error('build_adaptive_tv_RtR: inv_model缺少elem_data或jacobian_bkgnd.value');
end

% parameters for element-domain processing
avg_iters = 1;   % adjacency averaging iterations
lambda    = 1;   % self-weight in averaging

do_opening = (p < 0);  % not used normally
do_closing = (p >= 0);

% 1) edge gradients
L = prior_TV(imdl);                 % G x Ne
g_edge = abs(L * base_elem_data);   % G x 1

% 2) edge -> element aggregate
EInc = spones(abs(L));              % G x Ne incidence
edges_per_elem = full(max(1, sum(EInc,1)'));
g_elem = (EInc' * g_edge) ./ edges_per_elem; % Ne x 1

% 3) adjacency average on elements
A = spones(EInc' * EInc);
A = A - spdiags(diag(A),0,size(A,1),size(A,2));
deg = full(sum(A,2));
if avg_iters > 0
    for t = 1:avg_iters
        g_elem = (A * g_elem + lambda * g_elem) ./ max(1, (deg + lambda));
    end
end

% 4) morphology (optional)
if do_opening
    g_elem = element_graph_erosion(g_elem, A);
    g_elem = element_graph_dilation(g_elem, A);
end
if do_closing
    g_elem = element_graph_dilation(g_elem, A);
    g_elem = element_graph_erosion(g_elem, A);
end

% 5) element -> edge mapping (average of adjacent elements)
G = size(L,1);
g_edge_ref = zeros(G,1);
for r = 1:G
    elems = find(L(r,:) ~= 0);
    if isempty(elems)
        g_edge_ref(r) = 0;
    else
        g_edge_ref(r) = mean(g_elem(elems));
    end
end

% 6) produce weights and RtR
w = 1 ./ (1 + (g_edge_ref./max(tau,eps)).^p);
w = min(max(w, w_min), w_max);
W = spdiags(w, 0, numel(w), numel(w));
Lw = W * L;
RtR = Lw' * Lw;
end
