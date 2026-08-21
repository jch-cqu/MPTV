function img = apply_conductivity_pattern(img, sigma_ground, sigma_background, pattern_option,sigma_f)
% APPLY_CONDUCTIVITY_PATTERN  helper to set up ground/background strip pattern
% img = apply_conductivity_pattern(img, sigma_ground, sigma_background)
% img = apply_conductivity_pattern(..., pattern_option)
%
% pattern_option (optional):
%   3 - legacy grid/indices pattern (default)
%   4 - user-defined region mask pattern (uses interp_mesh)

if nargin < 4 || isempty(pattern_option)
    pattern_option = 3;
end

if pattern_option == 3
    % legacy behavior: set blocks and indices
    grid_blocks = [1 48; 385 432; 769 816; 1105 1152];
    for row = 1:size(grid_blocks, 1)
        img.elem_data(grid_blocks(row, 1):grid_blocks(row, 2)) = sigma_ground;
    end

    indices1 = arrayfun(@(x) x:x+1, 49:48:1058, 'UniformOutput', false);
    indices2 = arrayfun(@(x) x:x+1, 95:48:1104, 'UniformOutput', false);
    indices3 = arrayfun(@(x) x:x+1, 65:48:1074, 'UniformOutput', false);
    indices4 = arrayfun(@(x) x:x+1, 81:48:1090, 'UniformOutput', false);
    indices = [indices1, indices2, indices3, indices4];
    try
        img.elem_data([indices{:}]) = sigma_ground;
    catch
        % ignore if indices out of range for different meshes
    end

    % img.elem_data([1121:1136])=sigma_f;
    % img.elem_data([785:802])=sigma_f;%%9
      % img.elem_data([401:418])=sigma_f;%%16
     % img.elem_data([1105:1120])=sigma_f;%%1
     % img.elem_data([1121:1136])=sigma_f;%%2

%%
% 14支路
indices7 = arrayfun(@(x) x:x+1, 479:48:768, 'UniformOutput', false);
indices7 = cell2mat(indices7); % 将元胞数组转换为普通数组 
img.elem_data(indices7)  =0.5; % 使用转换后的数组进行索引 

%%
% 11支路
indices7 = arrayfun(@(x) x:x+1, 433:48:721, 'UniformOutput', false);
indices7 = cell2mat(indices7); % 将元胞数组转换为普通数组 
img.elem_data(indices7)  = sigma_f; % 使用转换后的数组进行索引 

%%
     %%%%%%段
    % img.elem_data([775:780])=sigma_f;%%8半
    % img.elem_data([791:796])=sigma_f;%%9半

    % img.elem_data([783:788])=sigma_f;%%横中
   %%竖中
    % indices6 = arrayfun(@(x) x:x+1, 737:48:833, 'UniformOutput', false);
    % indices6 = cell2mat(indices6); % 将元胞数组转换为普通数组 
    % img.elem_data(indices6)  = sigma_f; % 使用转换后的数组进行索引

elseif pattern_option == 2
    % Pattern taken from cc22.m: a custom tiled/strip layout using interp_mesh
    % Start from explicit background
    img.elem_data(:) = sigma_background;

    xym = interp_mesh(img, 3);
    x_xym = xym(:,1,:); y_xym = xym(:,2,:);

    ff = (x_xym>-1) & (x_xym<1) & (y_xym<-0.833) & (y_xym>-1);
    img.elem_data = img.elem_data + 4*mean(ff,3);
    ff = (x_xym>-1) & (x_xym<1) & (y_xym<0.0833) & (y_xym>-0.083);
    img.elem_data = img.elem_data + 4*mean(ff,3);
    ff = (x_xym>-1) & (x_xym<1) & (y_xym<1) & (y_xym>0.833);
    img.elem_data = img.elem_data + 4*mean(ff,3);
    ff = (x_xym>-1) & (x_xym<-0.83) & (y_xym<1) & (y_xym>-1);
    img.elem_data = img.elem_data + 4*mean(ff,3);
    ff = (x_xym>0.83) & (x_xym<1) & (y_xym<1) & (y_xym>-1);
    img.elem_data = img.elem_data + 4*mean(ff,3);
    ff = (x_xym>-0.0833) & (x_xym<0.083) & (y_xym<1) & (y_xym>-1);
    img.elem_data = img.elem_data + 4*mean(ff,3);

    % apply same clipping/mapping strategy as cc22: group and map to final values
    img.elem_data(img.elem_data == 4 | img.elem_data > 4) = 4;
    img.elem_data(img.elem_data > 2) = 4;
    img.elem_data(img.elem_data < 2) = 1;
    % Map back to sigma values: flagged (4) -> sigma_ground, others -> sigma_background
    img.elem_data(img.elem_data == 4 ) = sigma_ground;
    img.elem_data(img.elem_data == 1 ) = sigma_background;

    %%%number 7
    img.elem_data([553:576])=sigma_f;
    img.elem_data([601:624])=sigma_f;
    %%%%number 2
    % img.elem_data([1081:1104])=sigma_f;
    % img.elem_data([1129:1152])=sigma_f;
    % 
        % indices6 = arrayfun(@(x) x:x+1, 119:48:1056, 'UniformOutput', false);
        % indices6 = cell2mat(indices6);
        % img.elem_data(indices6)  = sigma_background;
   

elseif pattern_option == 4
    % User-defined region mask pattern using interp_mesh (your custom layout)
    % Start from explicit background
    img.elem_data(:) = sigma_background;

    xym = interp_mesh(img, 3);
    x_xym = xym(:,1,:); y_xym = xym(:,2,:);

    ff = (x_xym>-1) & (x_xym<-0.9236) & (y_xym<1) & (y_xym>-1);
    img.elem_data = img.elem_data + 4*mean(ff,3);
    ff = (x_xym>-0.5) & (x_xym<-0.416) & (y_xym<1) & (y_xym>-1);
    img.elem_data = img.elem_data + 2*mean(ff,3);
    ff = (x_xym>0) & (x_xym<0.0833) & (y_xym<1) & (y_xym>-1);
    img.elem_data = img.elem_data + 4*mean(ff,3);
    ff = (x_xym>0.5) & (x_xym<0.5833) & (y_xym<1) & (y_xym>-1);
    img.elem_data = img.elem_data + 4*mean(ff,3);
    ff = (x_xym>0.9166) & (x_xym<1) & (y_xym<1) & (y_xym>-1);
    img.elem_data = img.elem_data + 4*mean(ff,3);
    ff = (x_xym>-1) & (x_xym<1) & (y_xym>-1) & (y_xym<-0.9236);
    img.elem_data = img.elem_data + 4*mean(ff,3);
    ff = (x_xym>-1) & (x_xym<1) & (y_xym<-0.416) & (y_xym>-0.5);
    img.elem_data = img.elem_data + 4*mean(ff,3);
    ff = (x_xym>-1) & (x_xym<1) & (y_xym<0.0833) & (y_xym>0);
    img.elem_data = img.elem_data + 4*mean(ff,3);
    ff = (x_xym>-1) & (x_xym<1) & (y_xym<0.5833) & (y_xym>0.5);
    img.elem_data = img.elem_data + 4*mean(ff,3);
    ff = (x_xym>-1) & (x_xym<1) & (y_xym<1) & (y_xym>0.9166);
    img.elem_data = img.elem_data + 4*mean(ff,3);

    % Normalize / map to final conductivity values consistent with original logic
    img.elem_data(img.elem_data == 4 | img.elem_data > 4) = 4;
    img.elem_data(img.elem_data > 2) = 4;
    img.elem_data(img.elem_data < 2) = 1;
    img.elem_data(img.elem_data == 4 ) = sigma_ground;
    img.elem_data(img.elem_data == 1 ) = sigma_background;

    % %%%%%%%%%16
        % indices6 = arrayfun(@(x) x:x+1, 601:48:841, 'UniformOutput', false);
        % indices6 = cell2mat(indices6);
        % img.elem_data(indices6)  = sigma_f;
    %%%%%%%%%7
   %      indices6 = arrayfun(@(x) x:x+1, 937:48:1130, 'UniformOutput', false);
   %      indices6 = cell2mat(indices6);
   %      img.elem_data(indices6)  = sigma_f;
   % % %%%%%%%%28
   img.elem_data([291:302])=sigma_f;
   %%%%%%%%%20
   %  img.elem_data([590:600])=sigma_f;
   % %%%%%%%%21
   %  img.elem_data([601:613])=sigma_f;
    %%%%%10
    img.elem_data([867:878])=sigma_f;
     %%%%%%11
    % img.elem_data([879:889])=sigma_f;
    % %%%%%%%%%%1
    % % img.elem_data([1107:1116])=sigma_f;
    % %%%%%%%%%2
    % img.elem_data([1118:1129])=sigma_f;
    %%%%%%%%%19
    % img.elem_data([577:590])=sigma_f;
    %%%%%%%%%22
    % img.elem_data([615:624])=0.4;

else
    error('apply_conductivity_pattern: unknown pattern_option %d', pattern_option);
end

 
end
