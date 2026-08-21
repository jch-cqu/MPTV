
clear;
imb=  mk_common_model('d2s',12);
% 调整模型的几何尺寸
% 
f=[601 609 617 625 425 225 25 17 9 1 201 401];%33 12电极
    for i=1:12
    imb.fwd_model.electrode(i).nodes=f(i);
    end

bkgnd= 0.6;
bkgnd2= 0.12;
% Solve Homogeneous model（求解均匀介质模型）
img= mk_image(imb.fwd_model, bkgnd2);
% show_fem(img);
vh= fwd_solve( img );
% 
% 设置接地网模型电导率
% img.elem_data([1:32])=bkgnd * 2;
% 3×3  d
 img.elem_data([1:48])=bkgnd ;
 img.elem_data([1105:1152])=bkgnd ;
 img.elem_data([385:432])=bkgnd ;
 img.elem_data([769:816])=bkgnd ;

indices1 = arrayfun(@(x) x:x+1, 49:48:1058, 'UniformOutput', false);
indices2 = arrayfun(@(x) x:x+1, 95:48:1104, 'UniformOutput', false);
indices3 = arrayfun(@(x) x:x+1, 65:48:1074, 'UniformOutput', false);
indices4 = arrayfun(@(x) x:x+1, 81:48:1090, 'UniformOutput', false);
indices = [indices1, indices2, indices3, indices4];
img.elem_data([indices{:}]) = bkgnd ;
% %腐蚀位置处设立
 % img.elem_data([785:802])=.12;
  % img.elem_data([802:814])=.12;
  % img.elem_data([401:418])=.12;
 % img.elem_data([1105:1121])=.35;
img.elem_data([1121:1136])=.35;

% 
% indices5 = arrayfun(@(x) x:x+1, 465:48:848, 'UniformOutput', false);
% indices5 = cell2mat(indices5); % 将元胞数组转换为普通数组 
% img.elem_data(indices5)  = .35; % 使用转换后的数组进行索引 

%%第5支路去掉
% indices6 = arrayfun(@(x) x:x+1, 833:48:1120, 'UniformOutput', false);
% indices6 = cell2mat(indices6); % 将元胞数组转换为普通数组 
% img.elem_data(indices6)  = .35; % 使用转换后的数组进行索引 
% %%第12支路去掉
% indices6 = arrayfun(@(x) x:x+1, 449:48:737, 'UniformOutput', false);
% indices6 = cell2mat(indices6); % 将元胞数组转换为普通数组 
% img.elem_data(indices6)  = .32 % 使用转换后的数组进行索引
% %%第13支路去掉
% indices6 = arrayfun(@(x) x:x+1, 465:48:753, 'UniformOutput', false);
% indices6 = cell2mat(indices6); % 将元胞数组转换为普通数组 
% img.elem_data(indices6)  = .1 % 使用转换后的数组进行索引
% 第19支路去掉
% indices6 = arrayfun(@(x) x:x+1, 65:48:353, 'UniformOutput', false);
% indices6 = cell2mat(indices6); % 将元胞数组转换为普通数组 
% img.elem_data(indices6)  = .3; % 使用转换后的数组进行索引
%%第5支路左移
% indices7 = arrayfun(@(x) x:x+1, 841:48:1128, 'UniformOutput', false);
% indices7 = cell2mat(indices7); % 将元胞数组转换为普通数组 
% img.elem_data(indices7)  = .6; % 使用转换后的数组进行索引 
%%第12支路右移动
% indices7 = arrayfun(@(x) x:x+1, 457:48:792, 'UniformOutput', false);
% indices7 = cell2mat(indices7); % 将元胞数组转换为普通数组 
% img.elem_data(indices7)  = 0.5; % 使用转换后的数组进行索引 
%%第19支路右移动
% indices7 = arrayfun(@(x) x:x+1, 73:48:457, 'UniformOutput', false);
% indices7 = cell2mat(indices7); % 将元胞数组转换为普通数组 
% img.elem_data(indices7)  = .6; % 使用转换后的数组进行索引 

vi= fwd_solve( img );
show_fem(img);


inv2d= eidors_obj('inv_model', 'EIT inverse');
inv2d.reconst_type= 'absolute';
inv2d.jacobian_bkgnd.value= 1;
 
% This is not an inverse crime; inv_mdl != fwd_mdl
imb=  mk_common_model('i2s',12);
% 
f=[9313,9345,9377,9409,6208,3201,97,65,33,1,3105,6209];%33 12电极
    for i=1:12
    imb.fwd_model.electrode(i).nodes=f(i);
    end

inv2d.fwd_model= imb.fwd_model;


% % Guass-Newton solvers
% inv2d.solve=       @inv_solve_core;
% inv2d.hyperparameter.func = @prior_tikhonov;
% % inv2d.hyperparameter.noise_figure= 0.8;
% % inv2d.hyperparameter.tgt_elems= 1:4;
% inv2d.RtR_prior= @prior_laplace;
% inv2d.solve= @inv_solve_core;
% imgr(4)= inv_solve( inv2d, vi);
% % 
% inv2d.hyperparameter = rmfield(inv2d.hyperparameter,'func');
% % 
%     % Total variation using PDIPM
%   inv2d.RtR_prior=   @prior_laplace_old;
%     imgr(3)= inv_solve( inv2d, vi);
% 
% inv2d.hyperparameter.value = .012;
% inv2d.solve= @inv_solve_core;
% inv2d.R_prior= @prior_TV;
% inv2d.inv_solve_core.max_iterations= 10;
% % inv2d.inv_solve_core.term_tolerance= 1e-10;    
% inv2d.inv_solve_core.verbose = 3;
% imgr(5)= inv_solve( inv2d, vi);
%     imgr(5).calc_colours.window_range=2;
%     % 创建并保存图像
    % figure;
    % show_slices(imgr(5));

inv2d.solve=       @inv_solve_core;
 inv2d.hyperparameter.value = .0015;
    inv2d.RtR_prior=   @prior_tikhonov;
    imgr(1)= inv_solve( inv2d,  vi);

 % inv2d.hyperparameter.value = .001;
 %    inv2d.RtR_prior=   @prior_noser;
 %    imgr(2)= inv_solve( inv2d, vi);
  % inv2d.hyperparameter.value = .05;
  %   inv2d.RtR_prior=   @prior_laplace_old;
  %   imgr(3)= inv_solve( inv2d, vi);
        figure;
    show_slices(imgr(1));