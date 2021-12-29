%%
%该部分的功能为提取与保存数据
%文件初始化
clear;
clc;
close all;
load CD_run_params
output_path = './';
data = zeros(1536,2048);


for b = 1 : 2
file_pre = strcat( output_path, output_prefix, '_', num2str(b) );
AGC_values = load_AGC_block( file_pre, first_rg_line, ...
Nrg_lines_blk, b , UseMATfiles ); 
data(:,:,b) = load_DATA_block( file_pre, output_path, ...
Nrg_lines_blk, Nrg_cells, AGC_values, b, UseMATfiles );
end
data = double(data);
s_echo = [data(:,:,2);data(:,:,1)];

figure;
image(abs(s_echo));
title('原始图像/二维时域');
colormap(gray);
colorbar;
%%
%该部分的功能为定义一些参数
%方位向：Azimuth
%距离向：Range
Kr = -Kr;                       % 将调频率Kr改成负值
BW_range = 30.111e+06;          % 脉冲宽度
Vr = 7062;                      % 有效雷达速率
KA = 1733;                      % 方位调频率
f_doppler = -6900;                    % 多普勒中心频率
FA = PRF;                       % 方位向采样率
lamda = c/f0;                   % 波长
T_start = 6.5959e-03;           % 数据窗开始时间
Nr = round(Tr*Fr);              % 线性调频信号采样点数
Nrg = Nrg_cells;                % 距离线采样点数
Naz = Nrg_lines;          % 总共的距离线数
NFFT_rg = Nrg;            % 距离向FFT长度
NFFT_az = Naz;           % 方位向FFT长度
% 距离
% 距离时间轴
tr = T_start + ( -Nrg/2 : (Nrg/2-1) )/Fr;                       
% 距离频率轴
fr = ( -NFFT_rg/2 : NFFT_rg/2-1 )*( Fr/NFFT_rg );                 
% 方位
 % 方位时间轴
ta = ( -Naz/2: Naz/2-1 )/FA; 
% 方位频率轴
fa = f_doppler + ...
fftshift(-NFFT_az/2:NFFT_az/2-1)/NFFT_az*FA; 

%%  
% 距离压缩
% 首先先进行距离向上的傅里叶变换
S_range = fft(s_echo,NFFT_rg,2);  
% 生成距离向匹配滤波器，采用方式2
% 时域复制脉冲，末端补零后进行FFT，最后取复共轭
t_ref = ( -Nr/2 : (Nr/2-1) )/Fr;    % 用来生成距离时间轴
t_ref_mtx = ones(Naz,1)*t_ref;      % 矩阵形式
w_ref = kaiser(Nr,2.5);             % 距离向，构建Kaiser窗
w_ref = ones(Naz,1)*(w_ref.');      
% 构成矩阵形式，每一行都进行加窗操作

s_ref = w_ref.*exp((1j*pi*Kr).*((t_ref_mtx).^2)); 
%加了窗的复制（发射）脉冲。
s_ref = [s_ref,zeros(Naz,Nrg-Nr)];     
% 对复制脉冲，后端补零。
 
S_ref = fft(s_ref,NFFT_rg,2);         
% 复制脉冲的距离傅里叶变换。
H_range = conj(S_ref);                 % 距离向匹配滤波器。
S_range_c = S_range.*H_range;  % 乘以匹配滤波器。          
% 为了作比较，先变换回二维时域
s_rc = ifft(S_range_c,[],2);   

% 作图显示
figure;
imagesc(abs(s_rc));
title('距离压缩后/二维时域'); 
colormap(gray);
colorbar;

%%
% 该部分完成的功能为二次距离压缩（SRC）
% 距离压缩需要在频域完成
% 在二维时域进行数据搬移
s_rc = s_rc.*exp(-1j*2*pi*f_doppler.*(ta.'*ones(1,Nrg))); 
% 方位向傅里叶变换，到距离多普勒域
S_2df = fft(s_rc,NFFT_az,1);
% 距离向傅里叶变换，到二维频域
S_2df = fft(S_2df,Nrg,2);               

% 大斜视角下的徙动因子
D_fn_Vr = sqrt(1-lamda^2.*(fa.').^2./(4*Vr^2));         
K_src = 2*Vr^2*f0^3.*D_fn_Vr.^3./(c*R0*(fa.').^2);   
K_src_1 = 1./K_src;
% 二次距离压缩滤波器
H_src = exp(-1j*pi.*K_src_1*(fr.^2)); 
H_src = fftshift(H_src,2);   

S_2df_src = S_2df.*H_src;  
% 完成二次距离压缩（SRC），回到距离多普勒域。
S_rd = ifft(S_2df_src,[],2);  


% 作图显示
figure;
imagesc(abs(S_rd));
title('SRC后/距离多普勒域');
colormap(gray);
colorbar
%%
% 该部分的功能为完成距离徙动校正
% 距离徙动校正应该在距离多普勒域完成
% 而在上一部分中已经将其从二维频域搬移到距离多普勒域
% 每一个最近斜距（R0）都随着距离门的不同而改变。
tr_RCMC = T_start + ( -Nrg/2 : (Nrg/2-1) )/Fr;   % 用于RCMC和方位MF的距离时间轴。
R0_RCMC = (c/2).*tr_RCMC;   % 随距离线变化的R0，记为R0_RCMC，用来计算RCM和Ka。
delta_Rrd_fn = ((1-D_fn_Vr)./D_fn_Vr)*R0_RCMC;      % 大斜视角下的RCM

num_range = c/(2*Fr);   % 一个距离采样单元，对应的长度
delta_Rrd_fn_num = delta_Rrd_fn./num_range; % 每一个方位向频率，其RCM对应的距离采样单元数

R = 8;  % sinc插值核长度
S_rd_rcmc = zeros(NFFT_az,Nrg); % 用来存放RCMC后的值
for p = 1 : NFFT_az
    for q = 1 : Nrg   % 此时距离向的长度是 (Nrg-Nr+1)=Nrg        
        delta_Rrd_fn_p = delta_Rrd_fn_num(p,q);
        Rrd_fn_p = q + delta_Rrd_fn_p;
        
        Rrd_fn_p = rem(Rrd_fn_p,Nrg);  % 由于RCM的长度会超过Nrg，所以这样处理一下。
        
        Rrd_fn_p_int = ceil(Rrd_fn_p);        % ceil，向上取整。
        ii = ( Rrd_fn_p-(Rrd_fn_p_int-R/2):-1:Rrd_fn_p-(Rrd_fn_p_int+R/2-1)  );        
        rcmc_sinc = sinc(ii);
        rcmc_sinc = rcmc_sinc/sum(rcmc_sinc);   % 插值核的归一化
               
        % 由于S_rd只有整数点取值，且范围有限。因此插值中要考虑它的取值溢出边界问题。
        % 这里我采取循环移位的思想，用来解决取值溢出问题。
        if (Rrd_fn_p_int-R/2) > Nrg    % 全右溢
            variable = (Rrd_fn_p_int-R/2-Nrg:1:Rrd_fn_p_int+R/2-1-Nrg);
        else
            if (Rrd_fn_p_int+R/2-1) > Nrg    % 部分右溢
                variable_1 = (Rrd_fn_p_int-R/2:1:Nrg);
                variable_2 = (1:1:Rrd_fn_p_int+R/2-1-Nrg);
                variable = [variable_1,variable_2];
            else
                if (Rrd_fn_p_int+R/2-1) < 1    % 全左溢（不可能发生，但还是要考虑）
                    variable = (Rrd_fn_p_int-R/2+Nrg:1:Rrd_fn_p_int+R/2-1+Nrg);
                else
                    if (Rrd_fn_p_int-R/2) < 1       % 部分左溢
                        variable_1 = (Rrd_fn_p_int-R/2+Nrg:1:Nrg);
                        variable_2 = (1:1:Rrd_fn_p_int+R/2-1);
                        variable = [variable_1,variable_2];
                    else
                        variable = (Rrd_fn_p_int-R/2:1:Rrd_fn_p_int+R/2-1);
                    end                    
                end
            end
        end   
        rcmc_S_rd = S_rd(p,variable);
        S_rd_rcmc(p,q) = sum( rcmc_sinc.*rcmc_S_rd );
    end
end
% S_rd_rcmc 就是RCMC后的距离多普勒域频谱。

figure;
imagesc(abs(S_rd_rcmc));
title('RCMC后/距离多普勒域');
colormap(gray);
colorbar;
%%
%该部分完成的功能为方位压缩
fa_azimuth_MF = fa;         % 方位频率轴，采用和RCMC中所用的频率轴相同。
Haz = exp(1j*4*pi.*(D_fn_Vr*R0_RCMC).*f0./c);   % 改进的方位向MF

S_rd_c = S_rd_rcmc.*Haz;            % 乘以匹配滤波器
s_ac = ifft(S_rd_c,[],1);       	% 完成方位压缩，变到图像域。结束。

%%
% 全部的信号处理过程已经完成。接下来进行成像。
% 设置动态显示范围后，再显示图像。
sout = (s_ac);
sout = abs(sout)/max(max(abs(sout)));
G = 20*log10(sout+eps);                         % db显示
clim = [-55,0];                                           % 动态显示范围
figure;
imagesc(((0:Nrg-1)+first_rg_cell)/Fr*c/2+R0,((0:Naz-1)+first_rg_line)/FA*Vr,G,clim);
axis xy;
title('RADARSAT-1数据/RD算法成像')
xlabel('距离(m)')
ylabel('方位(m)')
colormap(gray);
colorbar;






