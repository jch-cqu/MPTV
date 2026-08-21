clear; clc; eidors_cache('clear_all');

%% ================== 算法参数（仅6个） ==================
hp_morph = 0.05;           % 整体正则化强度

opts = struct( ...
    'ridge',3e-6,...    % RtR对角项 2e-6
    'w_floor', 0.1, ...   % 权重下限 (骨架上的惩罚权重)0.3
    'w_ceil', 1.0, ...     % 权重上限 (背景的惩罚权重)
    'blend_alpha', 1, ...  % 
    'hthr_pct', 50, ...    % 高阈百分位 (marker)
    'lthr_pct', 30 ....    % 低阈百分位 (mask)
);

%% 1. 构建正向模
sigma_ground = 0.6;
sigma_background = 0.2;
sigma_f = 0.1;
pattern_option = 4;

% drive_electrodes = [601,613,625,325,25,13,1,301];%2*2
% drive_electrodes = [601 609 617 625 425 225 25 17 9 1 201 401];%3*3
drive_electrodes = [601 607 613 619 625 475 325 175 25 19 13 7 1 151 301 451];%4*4 
% drive_electrodes = 'auto';  % 使用EIDORS自带的16电极均匀分布

% fwd_mdl = mk_common_model('d2s', 8);%2*2
% fwd_mdl = mk_common_model('d2s', 12);
fwd_mdl = mk_common_model('d2s', 16);

if ~strcmp(drive_electrodes, 'auto')
    fwd_mdl.fwd_model = assign_electrode_nodes(fwd_mdl.fwd_model, drive_electrodes);
end

img_bg   = mk_image(fwd_mdl.fwd_model, sigma_background);
img_true = apply_conductivity_pattern(img_bg, sigma_ground, sigma_background, pattern_option,sigma_f);

vh = fwd_solve(img_bg);
vi = fwd_solve(img_true);
% ----- Add Gaussian noise (configurable)
% -----%%%%3x3双对称的seed顺序为70-0，70-50，75-20，75-200,70-100
add_noise_flag = false;      % set false to run noiseless
target_snr_db = 70;         % desired SNR in dB (power dB)
% noise_seed: numeric scalar for reproducible noise, [] or 'shuffle' for non-fixed
noise_seed = 200;             % set to [] or 'shuffle' for non-fixed noise
if add_noise_flag
  % set RNG according to requested seed and capture used seed for reproducibility
  if isempty(noise_seed) || (ischar(noise_seed) && strcmp(noise_seed,'shuffle'))
    rng('shuffle');
  else
    rng(noise_seed);
  end
  s = rng; used_seed = s.Seed;

  % add_noise expects SNR as amplitude (norm(signal)/norm(noise)).
  % Convert target dB (power definition) to amplitude ratio:
  SNR_amp = 10^(target_snr_db/20); % amplitude ratio for target dB
  vi_clean = vi; % keep reference
  vi = add_noise(SNR_amp, vi);
  noise = vi.meas - vi_clean.meas;
  actual_snr = 20*log10(norm(vi_clean.meas)/norm(noise));
  vi.noise_seed = used_seed;   % store used RNG seed for reproducibility
  vi.noise = noise;            % store noise vector for inspection if needed
  fprintf('Added Gaussian noise (seed=%d): target SNR=%gdB, actual SNR=%0.2gdB\n', used_seed, target_snr_db, actual_snr);
end
figure; show_fem(img_true); title('Ground truth');

%% 2. 构建反演模型
%2x2
% measurement_electrodes = [9313,9361,9409,4753,97,49,1,4657];
% inv_base = mk_common_model('i2s', 8);
%3x3
% measurement_electrodes = [9313 9345 9377 9409 6208 3201 97 65 33 1 3105 6209];
% inv_base = mk_common_model('i2s', 12);
%4x4 手动指定
measurement_electrodes = [9313,9337,9361,9385,9409,7081,4753,2425,97,73,49,25,1,2329,4657,6985];
% measurement_electrodes = 'auto';  % 使用EIDORS自带的16电极均匀分布
inv_base = mk_common_model('i2s', 16);

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
% 传递电极节点信息以启用边界保护
opts.electrode_nodes = measurement_electrodes;
imdl_morph.RtR_prior = @(m) build_path_regularization_RtR_with_morph_v2(m, opts);
if isfield(imdl_morph, 'R_prior'); imdl_morph = rmfield(imdl_morph, 'R_prior'); end

% 诊断
[RtR_full, diag] = build_path_regularization_RtR_with_morph_v2(imdl_morph, opts);
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

%% 5. 对比结果
figure('Name','Reconstruction Comparison');
subplot(1,2,1); show_slices(imgr_tv); title('Standard TV');
subplot(1,2,2); show_slices(imgr_morph); title('Morph Path-TV');


