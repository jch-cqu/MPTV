% for conductivity = 1e6:1e6:3e6    
    clear;
    % model=mphopen('"C:\Users\19022\Desktop\Untitled.mph"');
    model=mphopen('"C:\Users\19022\Desktop\论文资料\实验代码\comsol仿真\22(16电极)\2x2 （8电极）.mph"');
    mphgeom(model);
    %mphinterp%读入任意点的仿真数据


    % % 设置mat2的电导率
    model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'1e5[S/m]'});
    model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'1e7[S/m]'});
    % % 背景
    % terminal_selections = [67,91,115,139,177,170,163,156,132,108,84,60,12,20,28,36]; % 端点边界域的列表
    % ground_selections = [91,115,139,177,170,163,156,132,108,84,60,12,20,28,36,67]; % 地面边界域的列表
    %3x3(16电极)
    % terminal_selections = [96,153,184,244,295,288,281,274,237,177,146,89,36,43,50,57]; % 端点边界域的列表
    % ground_selections = [153,184,244,295,288,281,274,237,177,146,89,36,43,50,57,96]; % 地面边界域的列表
    % 3x3有土地激励源*****************
    % terminal_selections = [115,172,209,275,326,319,312,305,268,202,165,108,49,56,63,70]; % 端点边界域的列表
    % ground_selections = [172,209,275,326,319,312,305,268,202,165,108,49,56,63,70,115]; % 地面边界域的列表
    %8*8
    % terminal_selections = [304 408 435 541 603 599 595 591 537 429 361 309 224 228 232 236]; % 端点边界域的列表
    % ground_selections =   [408 435 541 603 599 595 591 537 429 361 309 224 228 232 236 304]; % 地面边界域的列表
    % 4*4改电极
    % terminal_selections = [29 109 181 260 338 319 223 334 330 254 177 101 25 33 105 37]; % 端点边界域的列表
    % ground_selections =   [109 181 260 338 319 223 334 330 254 177 101 25 33 105 37 29];
   % 2*3
   % terminal_selections = [85 107 164 209 255 248 241 234 202 153 100 78 28 35 42 49 ]; % 端点边界域的列表
   %  ground_selections =   [107 164 209 255 248 241 234 202 153 100 78 28 35 42 49 85];
   % 4*4模拟板子
   % terminal_selections = [73 184 287 390 452 446 440 434 427 324 221 118 7 27 43 59]; % 端点边界域的列表
   % ground_selections =   [184 287 390 452 446 440 434 427 324 221 118 7 27 43 59 73];
    % 33(12电极)
   % terminal_selections = [57 111 185 259 252 245 238 178 104 36 43 50]; % 端点边界域的列表
   % ground_selections =   [111 185 259 252 245 238 178 104 36 43 50 57];
   % 33(12电极在7x7网格种)
   % terminal_selections = [267 353 479 597 590 583 576 472 346 246 253 260]; % 端点边界域的列表
   % ground_selections =   [353 479 597 590 583 576 472 346 246 253 260 267];
    % 22(8电极)
   terminal_selections = [42 84 152 145 138 77 28 35]; % 端点边界域的列表
   ground_selections =   [84 152 145 138 77 28 35 42];
    % 初始化一个空矩阵来存储所有结果数据
    all_data = [];
    
    % 循环遍历每一对边界域
    for i = 1:length(ground_selections)
        % 设置地面边界域
        model.component('comp1').physics('ec').feature('gnd1').selection.set([ground_selections(i)]);
    
        % 设置端点边界域
        model.component('comp1').physics('ec').feature('term1').selection.set([terminal_selections(i)]);
    
        % 运行仿真
        model.sol('sol1').runAll;
        model.result.numerical('pev1').setResult;
    
        % 读取表格数据
        tbl = mphtable(model,'tbl1');
        data = tbl.data;
    
        % 将新数据按行添加到 all_data 矩阵中
        all_data = [all_data; data]; % 假设 data 是一个行向量
    
    end
    
        % all_data(:,1) = [];% 删除第一列***
          new_order = [1 4 7 6 8 5 2 3];%3x3
        % new_order = [6 8 10 12 16 15 14 13 11 9 7 5 1 2 3 4];%3x3
        % new_order = [5 7 9 11 16 15 14 13 12 10 8 6 1 2 3 4];%4x4模拟板子
        % new_order = [4 6 8 9 12 11 10 7 5 3 1 2];%3x3(12电极 （在7x7里面也一样）)
        % new_order = [5 8 10 12 16 15 14 13 11 9 7 6 1 2 3 4];%8*8
        % new_order = [3 4 8 7 11 12 16 15 14 13 9 10 6 5 1 2];%3x3有土地，改电极，褶形
        % new_order = [2 7 9 12 16 13 10 15 14 11 8 5 1 3 6 4 ];%4*4改电极
        % new_order = [6 8 10 12 16 15 14 13 11 9 7 5 1 2 3 4 ];%2*3
        new_data = all_data(:, new_order);
        new_data = abs(new_data);
        disp(['仿真已经执行完毕']);
    %-----------------------------------------------------------------------------------
    % clear;
    % clc;
   % disp(['正在运行腐蚀电导率为 ' num2str(conductivity) ' 的成像中，请等待……']);
    % [X1,X2,X3]=xlsread('C:\Users\19022\Desktop\9e6.xlsx');
    X1=new_data;
    X1=X1';
    vi.time([1])=NaN;
    vi.name='solved by fwd_solve_1st_order';
    vi.type='data';
    for i=1:8;
        for j=1:8;
            if j == 8;
                X4(i,j)=X1(1,i)-X1(j,i);
                % X4(i,j)=-X1(1,i)+X1(j,i);
            else 
            X4(i,j)=X1(j+1,i)-X1(j,i);
            end
        end
    end
    
    for u=1:8;
        i=1;
       switch u
            case 1;
                for h=1:5;
               vi.meas(1,(u-1)*5+h)=X4(u,2+h);
                end
                
            case 2;
                for h=1:5;
                 vi.meas(1,(u-1)*5+h)=X4(u,3+h);
                end
                 
            case 3;
                for h=1:5;
                    if h>1;   
                 vi.meas(1,(u-1)*5+h)=X4(u,3+h);
                    else
                 vi.meas(1,(u-1)*5+h)=X4(u,i);
                 i=i+1;
                    end
                end
                  
            case 4;
                for h=1:5;
                    if h>2;
                 vi.meas(1,(u-1)*5+h)=X4(u,3+h);
                    else
                 vi.meas(1,(u-1)*5+h)=X4(u,i);
                 i=i+1;
                    end
                end
                  
                 case 5;
                for h=1:5;
                    if h>3;
                 vi.meas(1,(u-1)*5+h)=X4(u,3+h);
                    else
                 vi.meas(1,(u-1)*5+h)=X4(u,i);
                 i=i+1;
                    end
                end
                  
                 case 6;
                for h=1:5;
                    if h>4;
                 vi.meas(1,(u-1)*5+h)=X4(u,3+h);
                    else
                 vi.meas(1,(u-1)*5+h)=X4(u,i);
                 i=i+1;
                    end
                end
                  
                 case 7;
                for h=1:5;
                    if h>5;
                 vi.meas(1,(u-1)*5+h)=X4(u,3+h);
                    else
                 vi.meas(1,(u-1)*5+h)=X4(u,i);
                 i=i+1;
                    end
                end
                  case 8;
                for h=2:6;
                 vi.meas(1,(u-1)*5+h-1)=X4(u,h);
                    end
                end
                  
        end
                    
    vi.meas=vi.meas';
    % for i=1:208
    for i=1:40
        vi.meas(i,1)=100*(vi.meas(i,1));
    end
      % vi.meas = round(vi.meas, 3);
     % 创建逆问题模型
    % 指定目标信噪比（SNR）以dB为单位
    % target_snr_db = 34;
    % % 计算信号的功率
    % signal_power = mean(vi.meas.^2);
    % % 根据SNR计算噪声功率
    % % SNR_db = 10 * log10(signal_power / noise_power)
    % % 因此，noise_power = signal_power / 10^(SNR_db / 10)
    % noise_power = signal_power / 10^(target_snr_db / 10);
    % % 计算噪声的标准差（幅度）
    % % 噪声功率是噪声方差，因此标准差是噪声方差的平方根
    % noise_std = sqrt(noise_power);
    % % 生成与信号相同长度的高斯噪声
    % noise = noise_std * randn(size(vi.meas));
    % % 将噪声添加到信号中
    % vi.meas = vi.meas + noise;

    