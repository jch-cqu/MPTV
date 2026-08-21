clear;
imb= mk_common_model('d2s',16);
  f=[601 607 613 619 625 475 325 175 25 19 13 7 1 151 301 451]; %16电极
    for i=1:16
    imb.fwd_model.electrode(i).nodes=f(i);
    end
e= size(imb.fwd_model.elems,1);
bkgnd= 0.3;

% Solve Homogeneous model（求解均匀介质模型）
img= mk_image(imb.fwd_model, bkgnd);
vh= fwd_solve( img );
% 设置接地网电导率
% 2×2 c
% 4×4 d
xym= interp_mesh( imb, 3);
x_xym= xym(:,1,:); y_xym= xym(:,2,:);
ff = (x_xym>-1) & (x_xym<-0.9236) & (y_xym<1) & (y_xym>-1);
img.elem_data= img.elem_data + 4*mean(ff,3);
ff = (x_xym>-0.5) & (x_xym<-0.416) & (y_xym<1) & (y_xym>-1);
img.elem_data= img.elem_data + 2*mean(ff,3);
ff = (x_xym>0) & (x_xym<0.0833) & (y_xym<1) & (y_xym>-1);
img.elem_data= img.elem_data + 4*mean(ff,3);
ff = (x_xym>0.5) & (x_xym<0.5833) & (y_xym<1) & (y_xym>-1);
img.elem_data= img.elem_data + 4*mean(ff,3);
x_xym= xym(:,1,:); y_xym= xym(:,2,:);
ff = (x_xym>0.9166) & (x_xym<1) & (y_xym<1) & (y_xym>-1);
img.elem_data= img.elem_data + 4*mean(ff,3);
ff = (x_xym>-1) & (x_xym<1) & (y_xym>-1) & (y_xym<-0.9236);
img.elem_data= img.elem_data + 4*mean(ff,3);
ff = (x_xym>-1) & (x_xym<1) & (y_xym<-0.416) & (y_xym>-0.5);
img.elem_data= img.elem_data + 4*mean(ff,3);
ff = (x_xym>-1) & (x_xym<1) & (y_xym<0.0833) & (y_xym>0);
img.elem_data= img.elem_data + 4*mean(ff,3);
ff = (x_xym>-1) & (x_xym<1) & (y_xym<0.5833) & (y_xym>0.5);
img.elem_data= img.elem_data + 4*mean(ff,3);
ff = (x_xym>-1) & (x_xym<1) & (y_xym<1) & (y_xym>0.9166);
img.elem_data= img.elem_data + 4*mean(ff,3);
img.elem_data(img.elem_data == 4 | img.elem_data > 4) = 4;
img.elem_data(img.elem_data > 2) = 4;
img.elem_data(img.elem_data < 2) = 1;
img.elem_data(img.elem_data == 4 ) = 0.5;
img.elem_data(img.elem_data == 1 ) = 0.3;

indices6 = arrayfun(@(x) x:x+1, 601:48:841, 'UniformOutput', false);
indices6 = cell2mat(indices6); % 将元胞数组转换为普通数组 
img.elem_data(indices6)  = .3; % 使用转换后的数组进行索引
% img.elem_data([291:302])=.3;
% img.elem_data([601:614])=.3;
% img.elem_data([867:878])=.3;
 % img.elem_data([1107:1116])=.3;
vi= fwd_solve( img );

show_fem(img);

inv2d= eidors_obj('inv_model', 'EIT inverse');
inv2d.reconst_type= 'absolute';
inv2d.jacobian_bkgnd.value= 1;
 
% This is not an inverse crime; inv_mdl != fwd_mdl
imb=  mk_common_model('i2s',16);

  f=[9313,9337,9361,9385,9409,7081,4753,2425,97,73,49,25,1,2329,4657,6985]; %16电极
    for i=1:16
    imb.fwd_model.electrode(i).nodes=f(i);
    end

inv2d.fwd_model= imb.fwd_model;
inv2d.solve=       @inv_solve_core;
inv2d.hyperparameter.func = @prior_tikhonov;
inv2d.hyperparameter.noise_figure= 0.8;
% inv2d.hyperparameter.tgt_elems= 1:4;
inv2d.RtR_prior= @prior_laplace;
inv2d.solve= @inv_solve_core;
imgr(4)= inv_solve( inv2d, vi);
% 
inv2d.hyperparameter = rmfield(inv2d.hyperparameter,'func');
% 
    % Total variation using PDIPM
inv2d.hyperparameter.value = .001;
inv2d.solve= @inv_solve_core;
inv2d.R_prior= @prior_TV;
inv2d.parameters.max_iterations= 12;
inv2d.parameters.term_tolerance= 1e-6;    
    imgr5= inv_solve( inv2d, vi);
    imgr5=rmfield(imgr5,'type'); 
    imgr5.type='image';
    imgr(5)=imgr5;
    % imgr(5).calc_colours.window_range=2;
    % 创建并保存图像
    figure;
    show_slices(imgr(5));
