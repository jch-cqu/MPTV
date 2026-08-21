clear; clc; eidors_cache('clear_all');

%% ================== 算法参数 ==================
hp_laplace = 0.01;           % Laplace正则化强度

%% 1. 构建正向模
sigma_ground = 0.6;
sigma_background = 0.2;
sigma_f = 0.4;
pattern_option = 2               ;

drive_electrodes = [601,613,625,325,25,13,1,301];%2*2
% drive_electrodes = [601 609 617 625 425 225 25 17 9 1 201 401];%3*3
% drive_electrodes = [601 607 613 619 625 475 325 175 25 19 13 7 1 151 301 451];%4*4

fwd_mdl = mk_common_model('d2s', 8);%2*2
% fwd_mdl = mk_common_model('d2s', 12);
% fwd_mdl = mk_common_model('d2s', 16);
fwd_mdl.fwd_model = assign_electrode_nodes(fwd_mdl.fwd_model, drive_electrodes);

img_bg   = mk_image(fwd_mdl.fwd_model, sigma_background);
img_true = apply_conductivity_pattern(img_bg, sigma_ground, sigma_background, pattern_option,sigma_f);

vh = fwd_solve(img_bg);
vi = fwd_solve(img_true);

figure; show_fem(img_true); title('Ground truth');

%% 2. 构建反演模型
%2x2
measurement_electrodes = [9313,9361,9409,4753,97,49,1,4657];
inv_base = mk_common_model('i2s', 8);
%3x3
% measurement_electrodes = [9313 9345 9377 9409 6208 3201 97 65 33 1 3105 6209];
% inv_base = mk_common_model('i2s', 12);
%4x4
% measurement_electrodes = [9313,9337,9361,9385,9409,7081,4753,2425,97,73,49,25,1,2329,4657,6985];
% inv_base = mk_common_model('i2s', 16);

inv_base.fwd_model = assign_electrode_nodes(inv_base.fwd_model, measurement_electrodes);

inv_template = eidors_obj('inv_model', 'Reconstruction');
inv_template.reconst_type = 'absolute';
inv_template.jacobian_bkgnd.value = 1; 
inv_template.fwd_model = inv_base.fwd_model;
inv_template.solve = @inv_solve_core;
inv_template.inv_solve_core.max_iterations = 10;
inv_template.inv_solve_core.term_tolerance = 1e-5;
inv_template.inv_solve_core.verbose = 1;

%% 3. Laplace 重建
imdl_laplace = inv_template;
imdl_laplace.hyperparameter.value = hp_laplace;
imdl_laplace.RtR_prior = @prior_laplace; 
if isfield(imdl_laplace, 'R_prior'); imdl_laplace = rmfield(imdl_laplace, 'R_prior'); end
imgr_laplace = inv_solve(imdl_laplace, vi);

%% 4. 结果展示
figure('Name','Laplace Reconstruction');
show_slices(imgr_laplace); 
title('Laplace');
colorbar;
