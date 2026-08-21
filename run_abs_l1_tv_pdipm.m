clear; clc; eidors_cache('clear_all');

%% ================== 算法参数 ==================
% Absolute L1-TV: L1 范数先验（TV/梯度稀疏），用 PDIPM 在绝对重建中求解
hp_tv = 0.01;              % 正则权重（建议从 1e-6 ~ 1e-3 扫描）
max_iter = 5;             % 最大迭代次数
beta = 1e-6;               % |x| 的平滑参数：abs(x) ≈ sqrt(x^2 + beta)

%% 1. 构建正向模
sigma_ground = 0.6;
sigma_background = 0.2;
sigma_f = 0.2;
pattern_option = 3;

% 2x2 / 3x3 / 4x4 可自行切换
% drive_electrodes = [601,613,625,325,25,13,1,301]; % 2*2
drive_electrodes = [601 609 617 625 425 225 25 17 9 1 201 401]; % 3*3
% drive_electrodes = [601 607 613 619 625 475 325 175 25 19 13 7 1 151 301 451]; % 4*4

% fwd_mdl = mk_common_model('d2s', 8);
fwd_mdl = mk_common_model('d2s', 12);
% fwd_mdl = mk_common_model('d2s', 16);
% If `drive_electrodes` is provided (non-empty) assign electrode nodes,
% otherwise keep the model's default electrode arrangement.
if exist('drive_electrodes','var') && ~isempty(drive_electrodes)
    fwd_mdl.fwd_model = assign_electrode_nodes(fwd_mdl.fwd_model, drive_electrodes);
else
    eidors_msg('Using model default electrode arrangement',2);
end

img_bg   = mk_image(fwd_mdl.fwd_model, sigma_background);
img_true = apply_conductivity_pattern(img_bg, sigma_ground, sigma_background, pattern_option, sigma_f);

vi = fwd_solve(img_true);

figure; show_fem(img_true); title('Ground truth');

%% 2. 构建反演模型（absolute + L1-TV(PDIPM)）
% 2x2
% measurement_electrodes = [9313,9361,9409,4753,97,49,1,4657];
% inv_base = mk_common_model('i2s', 8);
% 3x3
measurement_electrodes = [9313 9345 9377 9409 6208 3201 97 65 33 1 3105 6209];
inv_base = mk_common_model('i2s', 12);
% 4x4
% measurement_electrodes = [9313,9337,9361,9385,9409,7081,4753,2425,97,73,49,25,1,2329,4657,6985];
% inv_base = mk_common_model('i2s', 16);

% If `measurement_electrodes` provided, assign electrode nodes; otherwise
% use model defaults.
if exist('measurement_electrodes','var') && ~isempty(measurement_electrodes)
	inv_base.fwd_model = assign_electrode_nodes(inv_base.fwd_model, measurement_electrodes);
else
	eidors_msg('Using inverse model default electrode arrangement',2);
end

imdl = eidors_obj('inv_model', 'Absolute L1-TV (PDIPM) Reconstruction');
imdl.reconst_type = 'absolute';
imdl.fwd_model = inv_base.fwd_model;
imdl.jacobian_bkgnd.value = sigma_background;

% PDIPM 这里假定未知量与 fwd_model.elem_data 同维。
% 若模型包含 coarse2fine（粗到细映射），calc_R_prior/prior_TV 可能会把先验
% 映射到“粗参数空间”，从而导致 L*elem_data 维度不匹配。
if isfield(imdl.fwd_model, 'coarse2fine')
	imdl.fwd_model = rmfield(imdl.fwd_model, 'coarse2fine');
end

% PDIPM 绝对重建求解器（EIDORS 自带）
imdl.solve = @inv_solve_abs_pdipm;

% 数据项用 L2；先验项用 L1（TV）
imdl.inv_solve_abs_pdipm.norm_data  = 2;
imdl.inv_solve_abs_pdipm.norm_image = 1;
imdl.inv_solve_abs_pdipm.beta = beta;

imdl.hyperparameter.value = hp_tv;
imdl.parameters.max_iterations = max_iter;
imdl.parameters.min_change = 1e-6;

% TV 先验：R_prior = L（梯度/边算子），L1 发生在 |L*s| 上
% 使用数值型 TV 算子（避免缓存/参数化导致的维度漂移）
imdl.prior_use_fwd_not_rec = 1;
imdl.R_prior = prior_TV(imdl);
if isfield(imdl, 'RtR_prior'); imdl = rmfield(imdl, 'RtR_prior'); end

% 运行前做一次维度检查，便于定位问题
img_bkgnd_chk = calc_jacobian_bkgnd(imdl);
L_chk = calc_R_prior(imdl);
if size(L_chk,2) ~= size(img_bkgnd_chk.elem_data,1)
	error(['L1-TV(PDIPM) prior dimension mismatch: size(L,2)=%d but length(elem_data)=%d. ' ...
		   'Check coarse2fine/rec_model settings.'], size(L_chk,2), size(img_bkgnd_chk.elem_data,1));
end

%% 3. Absolute L1-TV 重建
imgr = inv_solve(imdl, vi);

%% 4. 结果展示
figure('Name','Absolute L1-TV (PDIPM) Reconstruction');
show_slices(imgr);
title('Absolute L1-TV (PDIPM)');
colorbar;
