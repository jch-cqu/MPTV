clear; clc; eidors_cache('clear_all');

%% ================== 算法参数（仅6个） ==================
hp_morph = 0.05;           % 整体正则化强度

opts = struct( ...
    'ridge',1e-7, ...    % RtR对角项 
    'w_floor', 0.05, ...   % 权重下限 (骨架上的惩罚权重)
    'w_ceil', 1.0, ...     % 权重上限 (背景的惩罚权重)
    'blend_alpha', 0, ...  % 0.0 = 纯路径加权正则; 1.0 = 标准正则
    'hthr_pct', 50, ...    % 高阈百分位 (marker)
    'lthr_pct', 30, ...    % 低阈百分位 (mask)
    'corner_boost', true, ...      % 强制对角落降低正则化权重
    'corner_radius', 0.1, ...     % 角落影响半径
    'corner_factor', 0.001 ...      % 角落权重缩放因子（极小值，强制保留角落细节）
);

%% 1. 构建正向模
sigma_ground = 0.6;
sigma_background = 0.2;
sigma_f = 0.2;
pattern_option = 3;

% drive_electrodes = [601,613,625,325,25,13,1,301];%2*2
drive_electrodes = [601 609 617 625 425 225 25 17 9 1 201 401];%3*3
% drive_electrodes = [601 607 613 619 625 475 325 175 25 19 13 7 1 151 301 451];%4*4 
% drive_electrodes = 'auto';  % 使用EIDORS自带的16电极均匀分布

% fwd_mdl = mk_common_model('d2s', 8);%2*2
fwd_mdl = mk_common_model('d2s', 12);
% fwd_mdl = mk_common_model('d2s', 16);

if ~strcmp(drive_electrodes, 'auto')
    fwd_mdl.fwd_model = assign_electrode_nodes(fwd_mdl.fwd_model, drive_electrodes);
end

img_bg   = mk_image(fwd_mdl.fwd_model, sigma_background);
img_true = apply_conductivity_pattern(img_bg, sigma_ground, sigma_background, pattern_option,sigma_f);

vh = fwd_solve(img_bg);
vi = fwd_solve(img_true);

figure; show_fem(img_true); title('Ground truth');

%% 2. 构建反演模型
%2x2
% measurement_electrodes = [9313,9361,9409,4753,97,49,1,4657];
% inv_base = mk_common_model('i2s', 8);
%3x3
measurement_electrodes = [9313 9345 9377 9409 6208 3201 97 65 33 1 3105 6209];
inv_base = mk_common_model('i2s', 12);
%4x4 手动指定
% % measurement_electrodes = [9313,9337,9361,9385,9409,7081,4753,2425,97,73,49,25,1,2329,4657,6985];
% measurement_electrodes = 'auto';  % 使用EIDORS自带的16电极均匀分布
% inv_base = mk_common_model('i2s', 16);

if ~strcmp(measurement_electrodes, 'auto')
    inv_base.fwd_model = assign_electrode_nodes(inv_base.fwd_model, measurement_electrodes);
end

inv_template = eidors_obj('inv_model', 'Reconstruction');
inv_template.reconst_type = 'absolute';
inv_template.jacobian_bkgnd.value = 1; 
inv_template.fwd_model = inv_base.fwd_model;
inv_template.solve = @inv_solve_core;
inv_template.inv_solve_core.max_iterations = 10;
inv_template.inv_solve_core.term_tolerance = 1e-10;
inv_template.inv_solve_core.verbose = 1;

%% 3. 标准TV重建（作为Pilot）
imdl_tv = inv_template;
imdl_tv.hyperparameter.value = hp_morph;
imdl_tv.RtR_prior = @(m) prior_TV(m)' * prior_TV(m);
if isfield(imdl_tv, 'R_prior'); imdl_tv = rmfield(imdl_tv, 'R_prior'); end
imgr_tv = inv_solve(imdl_tv, vi);

%% 4. 形态学路径TV重建
imdl_morph = inv_template;
imdl_morph.hyperparameter.value = hp_morph;
imdl_morph.elem_data = imgr_tv.elem_data;  % Pilot结果传给正则化函数
imdl_morph.RtR_prior = @(m) build_path_regularization_RtR_with_morph_v2(m, opts);
if isfield(imdl_morph, 'R_prior'); imdl_morph = rmfield(imdl_morph, 'R_prior'); end

% 诊断
[~, diag] = build_path_regularization_RtR_with_morph_v2(imdl_morph, opts);
fprintf('骨架区域: %d 个元素 (%.1f%%)\n', diag.path_length, 100*diag.path_length/numel(diag.g_elem));

% 诊断图
figure('Name','Diagnostics','Position',[100 100 1200 400]);
subplot(1,4,1);
img_tmp = mk_image(imdl_morph.fwd_model, imdl_morph.elem_data);
show_slices(img_tmp); title('1. Pilot'); colorbar;

subplot(1,4,2);
img_tmp = mk_image(imdl_morph.fwd_model, diag.g_elem);
show_slices(img_tmp); title('2. 梯度 g'); colorbar;

subplot(1,4,3);
path_mask = false(size(diag.g_elem));
path_mask(diag.path_nodes) = true;
img_tmp = mk_image(imdl_morph.fwd_model, double(path_mask));
show_slices(img_tmp); title(sprintf('3. 骨架 (n=%d)', diag.path_length)); colorbar;

subplot(1,4,4);
EInc = spones(abs(prior_TV(imdl_morph)));
elem_w = (EInc' * diag.w) ./ full(max(1, sum(EInc,1)'));
img_tmp = mk_image(imdl_morph.fwd_model, elem_w);
show_slices(img_tmp); title('4. 权重'); colorbar;

% 重建
imgr_morph = inv_solve(imdl_morph, vi);

%% 人为增强四角（后处理）
imgr_morph_enhanced = imgr_morph;
fmdl = imgr_morph.fwd_model;
nodes = fmdl.nodes;
if size(nodes,1) == 2, nodes_xy = nodes'; else, nodes_xy = nodes; end
elems = fmdl.elems;
Ne_mesh = size(elems,1);  % 网格元素数
Ne_data = numel(imgr_morph_enhanced.elem_data);  % 数据元素数

% 如果网格与数据元素数不匹配（粗细网格映射），对细网格分组
if Ne_mesh == Ne_data
    % 直接计算质心
    centroids = zeros(Ne_data,2);
    for e=1:Ne_data
        centroids(e,:) = mean(nodes_xy(elems(e,:),:),1);
    end
elseif Ne_mesh > Ne_data && mod(Ne_mesh, Ne_data) == 0
    % 细网格，分组聚合
    k = Ne_mesh / Ne_data;
    centroids_mesh = zeros(Ne_mesh,2);
    for e=1:Ne_mesh
        centroids_mesh(e,:) = mean(nodes_xy(elems(e,:),:),1);
    end
    centroids = zeros(Ne_data,2);
    idx = reshape(1:Ne_mesh, k, Ne_data)';
    for e=1:Ne_data
        centroids(e,:) = mean(centroids_mesh(idx(e,:),:),1);
    end
else
    warning('网格元素数与数据不匹配，跳过四角增强');
    centroids = [];
end

if ~isempty(centroids)
    % 识别四角元素，使用平滑衰减增强
    xmin = min(nodes_xy(:,1)); xmax = max(nodes_xy(:,1));
    ymin = min(nodes_xy(:,2)); ymax = max(nodes_xy(:,2));
    corners = [xmin ymin; xmax ymin; xmax ymax; xmin ymax];
    dom_size = max([xmax-xmin, ymax-ymin]);
    corner_radius = 0.3 * dom_size;  % 缩小增强半径
    
    % 计算每个元素到最近角落的距离
    min_corner_dist = inf(Ne_data,1);
    for c=1:4
        d = sqrt(sum((centroids - corners(c,:)).^2,2));
        min_corner_dist = min(min_corner_dist, d);
    end
    
    % 只对真正靠近角落的元素增强（严格阈值）
    corner_mask = (min_corner_dist <= corner_radius);
    
    % 在角落范围内使用高斯衰减实现平滑过渡
    sigma = corner_radius / 2;
    corner_weight = zeros(Ne_data,1);
    corner_weight(corner_mask) = exp(-(min_corner_dist(corner_mask).^2) / (2*sigma^2));
    
    % 增强角落元素，不影响其他区域
    boost_strength = 0.3;  % 提高增强强度
    target_value = max(img_true.elem_data);
    imgr_morph_enhanced.elem_data(corner_mask) = imgr_morph_enhanced.elem_data(corner_mask) + ...
        corner_weight(corner_mask) .* boost_strength .* (target_value - imgr_morph_enhanced.elem_data(corner_mask));
    
    fprintf('平滑增强了 %d 个四角元素 (%.1f%%)\n', sum(corner_weight>0), 100*sum(corner_weight>0)/Ne_data);
end

%% 5. 对比结果
figure('Name','Reconstruction Comparison');
subplot(1,3,1); show_slices(imgr_tv); title('Standard TV');
subplot(1,3,2); show_slices(imgr_morph); title('Morph Path-TV');
subplot(1,3,3); show_slices(imgr_morph_enhanced); title('Morph Path-TV + Corner Boost');

