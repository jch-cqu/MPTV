clear;clc;eidors_cache('clear_all');


% ============================================================================

% --- 主调参（放在脚本顶部，便于快速修改和网格搜索） ---
% 说明：把常做的调参暴露出来，数值会传入 build_path_regularization_RtR_with_morph。
% 优先调整 hp_morph（整体正则化强度），次调路径强度与阈值。
hp_morph    = 0.05;   % 推荐起点（你之前认为 0.025 效果好）
%若“没有路径”或 path_length 很小：先放宽 thr_factor（或降低 hthr_pct），然后放宽 morph_thresh（增加 mask 范围）
%thr_factor 与 morph_thresh 在源码中只在没有使用百分位策略（use_percentile==false）时才被用到；
thr_factor  = 0.30;    % marker / 强阈系数（决定哪些元素是候选路径）
morph_thresh = 0.20;   % mask 的低阈值系数（或用 lthr_pct 替代）
reduce_factor= 0.25;   % 路径衰减或降低因子（控制路径附近惩罚变小的幅度）

% 你仍可以通过 safe_opts 调整实现/稳定性选项（ridge, w_floor, blend_alpha 等）
% ================== 可调安全选项 (直接在脚本顶部修改以便试验) ==================
% 这些参数传入 build_path_regularization_RtR_with_morph 用于稳定 RtR 构造。
% 调参建议：先只改一个字段（例如 'blend_alpha' 或 'w_floor'），观察效果。
% 若“路径找到了但重建不稳定 / 数值警告”：先增 ridge，再考虑提高 w_floor。


 safe_opts = struct( ...
	'ridge', 1e-6, ...           % RtR 小对角项（提高w_floor后可用较小值）
	'w_floor', 0.25, ...          % 权重下限（提高以改善数值稳定性）
	'blend_alpha', 0.4, ...      % blend= a*(L'L) + (1-a)*(Lw'Lw)，提高以增加稳定性
	'inner_iters', 1, ...         % IRLS 内迭代次数
	'alpha', 0.8, ...             % IRLS 衰减系数
	'use_percentile', true, ...   % 使用百分位阈值
	'hthr_pct', 70, ...           % 高阈百分位 - 降低以包含更多区域
	'lthr_pct', 60, ...           % 低阈百分位 - 降低以扩大mask范围
	'w_ceil', 1.0, ...            % 权重上限（降低以减小权重范围）
	'blend_uniform', 0.1, ...    % 混合到均值以避免孤立极端权重
	'combine_mode', 'blend', ...  % 'plain' or 'blend'
	'normalize_mode', 'trace', ...% 归一化模式 ('none'|'trace'|'fro')
	'use_midline_band', false, ... % 优先在中轴带中挑端点
	'midline_band', 0.12, ...     % 中轴带宽度（模型 Y 轴绝对距离）
	'max_path_frac', 1, ...       % 路径长度上限（元素比例），经验0.08，无限为1
	'max_path_abs', inf ...       % 路径长度上限（元素绝对数），经验140，无限为inf
	);

%% 1. 构建正向模型
sigma_ground = 0.6;
sigma_background = 0.2;
sigma_f = 0.3;


% %%%%%%%%%%%%%%%%%%%%
pattern_option = 4;

% drive_electrodes = [601,613,625,325,25,13,1,301];%2*2
  % drive_electrodes = [601 609 617 625 425 225 25 17 9 1 201 401];%3*3
drive_electrodes = [601 607 613 619 625 475 325 175 25 19 13 7 1 151 301 451];%4*4

% % fwd_mdl = mk_common_model('d2s', 8);%2*2
% fwd_mdl = mk_common_model('d2s', 12);
fwd_mdl = mk_common_model('d2s', 16);


fwd_mdl.fwd_model = assign_electrode_nodes(fwd_mdl.fwd_model, drive_electrodes);

img_bg   = mk_image(fwd_mdl.fwd_model, sigma_background);
img_true = apply_conductivity_pattern(img_bg, sigma_ground, sigma_background, pattern_option,sigma_f);

vh = fwd_solve(img_bg);
vi = fwd_solve(img_true);

figure;
show_fem(img_true); title('Ground truth conductivity');


%% 2. 构建反演模型
%2x2
% measurement_electrodes = [9313,9361,9409,4753,97,49,1,4657];
% inv_base = mk_common_model('i2s', 8);
%3x3
% measurement_electrodes = [9313 9345 9377 9409 6208 3201 97 65 33 1 3105 6209];
% inv_base = mk_common_model('i2s', 12);
%4x4
measurement_electrodes = [9313,9337,9361,9385,9409,7081,4753,2425,97,73,49,25,1,2329,4657,6985];
inv_base = mk_common_model('i2s', 16);



inv_base.fwd_model = assign_electrode_nodes(inv_base.fwd_model, measurement_electrodes);

inv_template = eidors_obj('inv_model', 'Normalized absolute reconstruction');
inv_template.reconst_type = 'absolute';
inv_template.jacobian_bkgnd.value = 1;
inv_template.fwd_model = inv_base.fwd_model;
inv_template.solve = @inv_solve_core;
inv_template.inv_solve_core.max_iterations = 10;
inv_template.inv_solve_core.term_tolerance = 1e-9;
inv_template.inv_solve_core.verbose = 1;

% 先做一次标准TV重建作为 Pilot（避免使用常数背景导致 g=0 而退化为标准TV）
imdl_tv = inv_template;
imdl_tv.hyperparameter = struct('value', 0.05);
imdl_tv.RtR_prior = @(m) prior_TV(m)' * prior_TV(m);
if isfield(imdl_tv, 'R_prior'); imdl_tv = rmfield(imdl_tv, 'R_prior'); end
imgr_tv = inv_solve(imdl_tv, vi);

%  形态学重建增强的 Path 正则化（用 pilot 的 elem_data 引导路径检测）

% 使用当前的复杂 build_path_regularization_RtR_with_morph 构造器
imdl_morph = inv_template;
% 使用顶部暴露的主调参
imdl_morph.hyperparameter = struct('value', hp_morph); % 关键：保持与 TV 相同的超参数以公平比较
imdl_morph.elem_data = imgr_tv.elem_data; % 关键：提供非平凡的元素分布用于路径检测

% 自动按网格计算并覆盖有效路径上限 (非破坏式) —— 打印诊断以便调试
nelems = size(imdl_morph.fwd_model.elems,1);
user_abs = safe_opts.max_path_abs;
user_frac = safe_opts.max_path_frac;
if isfinite(user_abs)
	eff_max = min(user_abs, ceil(user_frac * nelems));
else
	eff_max = max(8, ceil(user_frac * nelems));
end
fprintf('Auto max_path_abs: mesh elems=%d, max_path_frac=%.3g -> eff_max=%d (user_abs=%s)\n', nelems, user_frac, eff_max, mat2str(user_abs));
% 覆盖 safe_opts 以传入 builder（非破坏性）
safe_opts.max_path_abs = eff_max;

% 使用顶部定义的 safe_opts（已覆盖 max_path_abs）
imdl_morph.RtR_prior = @(inv_model) build_path_regularization_RtR_with_morph(inv_model, thr_factor, reduce_factor, morph_thresh, safe_opts);
if isfield(imdl_morph, 'R_prior'); imdl_morph = rmfield(imdl_morph, 'R_prior'); end

% 先构造一次 RtR 并显示诊断图，便于在真正运行 inv_solve 前检查
try
	[RtR_test, diagnostics] = build_path_regularization_RtR_with_morph(imdl_morph, thr_factor, reduce_factor, morph_thresh, safe_opts);
	fprintf('Diagnostics: path_len=%d, mask2=%d\n', diagnostics.path_length, sum(diagnostics.mask2));
	
	% ===== 完整诊断对比图（无网格线） =====
	figure('Name','Full Diagnostics','Position',[100 100 1400 800]);
	
	% 1. Pilot 重建结果（输入数据）
	subplot(2,4,1);
	img_tmp = mk_image(imdl_morph.fwd_model, imdl_morph.elem_data);
	img_tmp.fwd_model.mdl_slice_mapper.npx = 128; img_tmp.fwd_model.mdl_slice_mapper.npy = 128;
	show_slices(img_tmp); title('1. Pilot σ (输入)'); colorbar;
	
	% 2. g_elem（梯度幅值场）
	subplot(2,4,2);
	img_tmp = mk_image(imdl_morph.fwd_model, diagnostics.g_elem);
	img_tmp.fwd_model.mdl_slice_mapper.npx = 128; img_tmp.fwd_model.mdl_slice_mapper.npy = 128;
	show_slices(img_tmp); title('2. g=|L\sigma| (梯度)'); colorbar;
	
	% 3. recon（形态学重建结果 = 骨架区域）
	subplot(2,4,3);
	img_tmp = mk_image(imdl_morph.fwd_model, double(diagnostics.recon));
	img_tmp.fwd_model.mdl_slice_mapper.npx = 128; img_tmp.fwd_model.mdl_slice_mapper.npy = 128;
	show_slices(img_tmp); title('3. recon (形态学重建)'); colorbar;
	
	% 4. mask2（= recon）
	subplot(2,4,4);
	img_tmp = mk_image(imdl_morph.fwd_model, double(diagnostics.mask2));
	img_tmp.fwd_model.mdl_slice_mapper.npx = 128; img_tmp.fwd_model.mdl_slice_mapper.npy = 128;
	show_slices(img_tmp); title(sprintf('4. mask2 (hthr=%.3f)', diagnostics.hthr)); colorbar;
	
	% 5. path_nodes（最大连通分量 = 骨架）
	path_mask = false(size(diagnostics.g_elem));
	if ~isempty(diagnostics.path_nodes)
		path_mask(diagnostics.path_nodes) = true;
	end
	subplot(2,4,5);
	img_tmp = mk_image(imdl_morph.fwd_model, double(path_mask));
	img_tmp.fwd_model.mdl_slice_mapper.npx = 128; img_tmp.fwd_model.mdl_slice_mapper.npy = 128;
	show_slices(img_tmp); title(sprintf('5. 接地网 (n=%d)', diagnostics.path_length)); colorbar;
	
	% 6. 元素权重
	EInc = spones(abs(prior_TV(imdl_morph)));
	edges_per_elem = full(max(1, sum(EInc,1)'));
	elem_weight = (EInc' * diagnostics.w) ./ edges_per_elem;
	subplot(2,4,6);
	img_tmp = mk_image(imdl_morph.fwd_model, elem_weight);
	img_tmp.fwd_model.mdl_slice_mapper.npx = 128; img_tmp.fwd_model.mdl_slice_mapper.npy = 128;
	show_slices(img_tmp); title('6. elem weight'); colorbar;
	
	% 7. path_scale（路径缩放因子）
	if isfield(diagnostics, 'path_scale')
		path_scale_elem = (EInc' * diagnostics.path_scale) ./ edges_per_elem;
		subplot(2,4,7);
		img_tmp = mk_image(imdl_morph.fwd_model, path_scale_elem);
		img_tmp.fwd_model.mdl_slice_mapper.npx = 128; img_tmp.fwd_model.mdl_slice_mapper.npy = 128;
		show_slices(img_tmp); title('7. path\_scale'); colorbar;
	end
	
	% 8. 说明文字
	subplot(2,4,8);
	axis off;
	text(0.05, 0.95, '论文流程:', 'FontSize', 11, 'FontWeight', 'bold');
	text(0.05, 0.80, '1. g=|L\sigma| 梯度幅值', 'FontSize', 9);
	text(0.05, 0.65, '2. 形态学开闭滤波', 'FontSize', 9);
	text(0.05, 0.50, '3. 双阈值+形态学重建', 'FontSize', 9);
	text(0.05, 0.35, '4. 最大连通分量=骨架', 'FontSize', 9);
	text(0.05, 0.20, sprintf('hthr_pct=%d%%, lthr_pct=%d%%', safe_opts.hthr_pct, safe_opts.lthr_pct), 'FontSize', 9);
	text(0.05, 0.05, sprintf('骨架=%d/%d (%.1f%%)', diagnostics.path_length, numel(diagnostics.g_elem), 100*diagnostics.path_length/numel(diagnostics.g_elem)), 'FontSize', 9, 'Color', 'b');
	
catch ME
	warning('build_path_regularization_RtR_with_morph:diag','%s', ME.message);
	disp(getReport(ME));
end


% 现在执行用复杂先验的重建
imgr_morph = inv_solve(imdl_morph, vi);

figure('Name','Morph Path-TV vs Standard TV');
subplot(1,2,1); show_slices(imgr_tv); title('Standard TV (pilot)');
subplot(1,2,2); show_slices(imgr_morph); title('Morphological Path-TV');
