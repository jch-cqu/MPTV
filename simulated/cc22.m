
clear;
% Compare 2D algorithms
% $Id: tutorial120a.m 3273 2012-06-30 18:00:35Z aadler $
imb= mk_common_model('d2s',8);

  f=[601,613,625,325,25,13,1,301];%33 8电极
    for i=1:8
    imb.fwd_model.electrode(i).nodes=f(i);
    end


e= size(imb.fwd_model.elems,1);
bkgnd= 0.3;
% Solve Homogeneous model（求解均匀介质模型）
img= mk_image(imb.fwd_model, bkgnd);
vh= fwd_solve( img );
% 设置接地网电导率
% 2×2 
xym= interp_mesh( imb, 3);
x_xym= xym(:,1,:); y_xym= xym(:,2,:);
ff = (x_xym>-1) & (x_xym<1) & (y_xym<-0.833) & (y_xym>-1);
img.elem_data= img.elem_data + 4*mean(ff,3);
x_xym= xym(:,1,:); y_xym= xym(:,2,:);
ff = (x_xym>-1) & (x_xym<1) & (y_xym<0.0833) & (y_xym>-0.083);
img.elem_data= img.elem_data + 4*mean(ff,3);
x_xym= xym(:,1,:); y_xym= xym(:,2,:);
ff = (x_xym>-1) & (x_xym<1) & (y_xym<1) & (y_xym>0.833);
img.elem_data= img.elem_data + 4*mean(ff,3);
x_xym= xym(:,1,:); y_xym= xym(:,2,:);
ff = (x_xym>-1) & (x_xym<-0.83) & (y_xym<1) & (y_xym>-1);
img.elem_data= img.elem_data + 4*mean(ff,3);
ff = (x_xym>0.83) & (x_xym<1) & (y_xym<1) & (y_xym>-1);
img.elem_data= img.elem_data + 4*mean(ff,3);
ff = (x_xym>-0.0833) & (x_xym<0.083) & (y_xym<1) & (y_xym>-1);
img.elem_data= img.elem_data + 4*mean(ff,3);
img.elem_data(img.elem_data == 4 | img.elem_data > 4) = 4;
img.elem_data(img.elem_data > 2) = 2;
img.elem_data(img.elem_data < 2) = 1;
img.elem_data(img.elem_data == 2) = 0.7;
img.elem_data(img.elem_data == 1) = 0.2;

indices5 = arrayfun(@(x) x:x+3, 119:48:1056, 'UniformOutput', false);
indices6 = arrayfun(@(x) x:x+3, 551:48:1056, 'UniformOutput', false);
indices5 = cell2mat(indices5); % 将元胞数组转换为普通数组 
img.elem_data(indices5)  = .42; % 使用转换后的数组进行索引 
indices6 = cell2mat(indices6); % 将元胞数组转换为普通数组 
img.elem_data(indices6)  = .3; % 使用转换后的数组进行索引 
% img.elem_data([551:572])=.3;img.elem_data([599:620])=.3;
% img.elem_data([5:24])=.3;img.elem_data([53:72])=.3;
% img.elem_data([1109:1128])=.35;img.elem_data([1061:1080])=.32;
vi= fwd_solve( img );

% 指定目标信噪比（SNR）以dB为单位
target_snr_db = 34;
% 计算信号的功率
signal_power = mean(vi.meas.^2);
% 根据SNR计算噪声功率
% SNR_db = 10 * log10(signal_power / noise_power)
% 因此，noise_power = signal_power / 10^(SNR_db / 10)
noise_power = signal_power / 10^(target_snr_db / 10);
% 计算噪声的标准差（幅度）
% 噪声功率是噪声方差，因此标准差是噪声方差的平方根
noise_std = sqrt(noise_power);
% 生成与信号相同长度的高斯噪声
noise = noise_std * randn(size(vi.meas));
% 将噪声添加到信号中
vi.meas = vi.meas + noise;

show_fem(img);

inv2d= eidors_obj('inv_model', 'EIT inverse');
inv2d.reconst_type= 'absolute';
inv2d.jacobian_bkgnd.value= 1;
 
% This is not an inverse crime; inv_mdl != fwd_mdl
imb=  mk_common_model('i2s',8);

 f=[9313,9361,9409,4753,97,49,1,4657];%33 8电极
    for i=1:8
    imb.fwd_model.electrode(i).nodes=f(i);
    end

inv2d.fwd_model= imb.fwd_model;
 
% % Guass-Newton solvers
inv2d.solve=       @inv_solve_core;

% inv2d.hyperparameter.func = @prior_tikhonov;
% inv2d.hyperparameter.noise_figure= 0.8;
% % inv2d.hyperparameter.tgt_elems= 1:4;
% inv2d.RtR_prior= @prior_laplace;
% inv2d.solve= @inv_solve_core;
% imgr(4)= inv_solve( inv2d, vi);
% % 
% inv2d.hyperparameter = rmfield(inv2d.hyperparameter,'func');
% % 
%     % Total variation using PDIPM
% inv2d.hyperparameter.value = .001;
% inv2d.solve= @inv_solve_core;
% inv2d.R_prior= @prior_TV;
% inv2d.inv_solve_core.max_iterations= 10;
% % inv2d.parameters.term_tolerance= 1e-4;    
%     imgr5= inv_solve( inv2d, vi);
%     imgr5=rmfield(imgr5,'type'); 
%     imgr5.type='image';
%     imgr(5)=imgr5;
%     % imgr(5).calc_colours.window_range=2;
%     % 创建并保存图像
%     figure;
%     show_slices(imgr(5));

   inv2d.hyperparameter.value = .005;
    inv2d.RtR_prior=   @prior_tikhonov;
    imgr(1)= inv_solve( inv2d,  vi);
 % 
 % inv2d.hyperparameter.value = .001;
 %    inv2d.RtR_prior=   @prior_noser;
 %    imgr(2)= inv_solve( inv2d, vi);
  % inv2d.hyperparameter.value = .005;
  %   inv2d.RtR_prior=   @prior_laplace;
  %   imgr(3)= inv_solve( inv2d, vi);
        figure;
    show_slices(imgr(1));