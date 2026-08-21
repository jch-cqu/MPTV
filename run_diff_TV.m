

clear; clc; eidors_cache('clear_all');

%% 参数
% No noise: comparisons use clean simulated data

%% 构建正向模（使用你的参数：d2s, 8）
sigma_ground = 0.6;
sigma_background = 0.2;
sigma_f = 0.2;
pattern_option =4;

% drive_electrodes = [601,613,625,325,25,13,1,301];%2*2
% drive_electrodes = [601 609 617 625 425 225 25 17 9 1 201 401];%3*3
drive_electrodes = [601 607 613 619 625 475 325 175 25 19 13 7 1 151 301 451];%4*4

% fwd_mdl = mk_common_model('d2s', 8);%2*2
% fwd_mdl = mk_common_model('d2s', 12);
fwd_mdl = mk_common_model('d2s', 16);
fwd_mdl.fwd_model = assign_electrode_nodes(fwd_mdl.fwd_model, drive_electrodes);

img_bg   = mk_image(fwd_mdl.fwd_model, sigma_background);
img_true = apply_conductivity_pattern(img_bg, sigma_ground, sigma_background, pattern_option, sigma_f);

vh = fwd_solve(img_bg);
vi = fwd_solve(img_true);
figure;
show_fem(img_true); title('Ground truth conductivity');
% use clean data only
vi_n = vi; 

%% 创建逆问题模板（参考 run_standard_tv 的 inv_template 风格）
% Use I2S model for inverse problem (measurement electrode layout)
%2x2
% measurement_electrodes = [9313,9361,9409,4753,97,49,1,4657];
% inv_base = mk_common_model('i2s', 8);
%3x3
% measurement_electrodes = [9313 9345 9377 9409 6208 3201 97 65 33 1 3105 6209];
% inv_base = mk_common_model('i2s', 12);
% % 4x4
measurement_electrodes = [9313,9337,9361,9385,9409,7081,4753,2425,97,73,49,25,1,2329,4657,6985];
inv_base = mk_common_model('i2s', 16);
inv_base.fwd_model = assign_electrode_nodes(inv_base.fwd_model, measurement_electrodes);

% Guass-Newton solvers
inv2d = eidors_obj('inv_model', 'EIT inverse GN');
inv2d.reconst_type = 'difference';
inv2d.jacobian_bkgnd.value = 1;
inv2d.solve =       @inv_solve_diff_GN_one_step;
inv2d.fwd_model = inv_base.fwd_model;
% %% Tikhonov prior
% inv2d.hyperparameter.value = .001;
% inv2d.RtR_prior=   @prior_tikhonov;
% imgr(1)= inv_solve( inv2d, vh, vi);
% imgn(1)= inv_solve( inv2d, vh, vi_n);
% figure;
% show_slices(imgn(1));
%% NOSER prior
% inv2d.hyperparameter.value = .003;
% inv2d.RtR_prior=   @prior_noser;
% imgr(2)= inv_solve( inv2d, vh, vi);
% imgn(2)= inv_solve( inv2d, vh, vi_n);
% figure;
% show_slices(imgn(2));
%% Laplace image prior
% inv2d.hyperparameter.value = .001;
% inv2d.RtR_prior=   @prior_laplace;
% imgr(3)= inv_solve( inv2d, vh, vi);
% imgn(3)= inv_solve( inv2d, vh, vi_n);
% figure;
% show_slices(imgn(3));

%% TV 监视窗口（保存中间迭代）
invtv= eidors_obj('inv_model', 'EIT inverse TV monitor');
invtv.reconst_type= 'difference';
invtv.jacobian_bkgnd.value= 1;
invtv.hyperparameter.value = .0001;%0.001在2x2和3x3
invtv.inv_solve_TV_pdipm.alpha1 = 0.05;%0.15在2x2和3x3
invtv.solve=       @inv_solve_TV_pdipm;
invtv.R_prior=     @prior_TV;
invtv.parameters.term_tolerance= 1e-6;
invtv.parameters.keep_iterations= 1; % keep iterations

invtv.fwd_model= inv_base.fwd_model;
invtv.parameters.max_iterations= 12;
imgtv= inv_solve( invtv, vh, vi);

% 显示 TV 结果（不保存图片）
figure('Name','TV monitor'); clf;
show_slices(imgtv);
title('TV monitor (kept iterations)');

% plot slice profiles across kept iterations if available
% try
%   imgs = calc_slices(imgtv);
%   idx = [1,3,6,10,20];
%   figure('Name','TV slice profiles'); clf;
%   subplot(2,1,1);
%   plot(squeeze(imgs(:,32,idx)));
%   axis([1 size(imgs,1) -0.04 0.12]); set(gca,'XTickLabel',[]);
%   legend(arrayfun(@num2str, idx, 'UniformOutput', false));
% 
%   subplot(2,1,2);
%   plot(squeeze(imgs(32,:,idx)));
%   axis([1 size(imgs,2) -0.04 0.12]); set(gca,'XTickLabel',[]);
% catch
%   % if intermediate iterations not available, just skip plots
% end
% %% 