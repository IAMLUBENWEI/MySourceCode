%%
%�ò��ֵĹ���Ϊ��ȡ�뱣������
%�ļ���ʼ��
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
title('ԭʼͼ��/��άʱ��');
colormap(gray);
colorbar;
%%
%�ò��ֵĹ���Ϊ����һЩ����
%��λ��Azimuth
%������Range
Kr = -Kr;                       % ����Ƶ��Kr�ĳɸ�ֵ
BW_range = 30.111e+06;          % ������
Vr = 7062;                      % ��Ч�״�����
KA = 1733;                      % ��λ��Ƶ��
f_doppler = -6900;                    % ����������Ƶ��
FA = PRF;                       % ��λ�������
lamda = c/f0;                   % ����
T_start = 6.5959e-03;           % ���ݴ���ʼʱ��
Nr = round(Tr*Fr);              % ���Ե�Ƶ�źŲ�������
Nrg = Nrg_cells;                % �����߲�������
Naz = Nrg_lines;          % �ܹ��ľ�������
NFFT_rg = Nrg;            % ������FFT����
NFFT_az = Naz;           % ��λ��FFT����
% ����
% ����ʱ����
tr = T_start + ( -Nrg/2 : (Nrg/2-1) )/Fr;                       
% ����Ƶ����
fr = ( -NFFT_rg/2 : NFFT_rg/2-1 )*( Fr/NFFT_rg );                 
% ��λ
 % ��λʱ����
ta = ( -Naz/2: Naz/2-1 )/FA; 
% ��λƵ����
fa = f_doppler + ...
fftshift(-NFFT_az/2:NFFT_az/2-1)/NFFT_az*FA; 

%%  
% ����ѹ��
% �����Ƚ��о������ϵĸ���Ҷ�任
S_range = fft(s_echo,NFFT_rg,2);  
% ���ɾ�����ƥ���˲��������÷�ʽ2
% ʱ�������壬ĩ�˲�������FFT�����ȡ������
t_ref = ( -Nr/2 : (Nr/2-1) )/Fr;    % �������ɾ���ʱ����
t_ref_mtx = ones(Naz,1)*t_ref;      % ������ʽ
w_ref = kaiser(Nr,2.5);             % �����򣬹���Kaiser��
w_ref = ones(Naz,1)*(w_ref.');      
% ���ɾ�����ʽ��ÿһ�ж����мӴ�����

s_ref = w_ref.*exp((1j*pi*Kr).*((t_ref_mtx).^2)); 
%���˴��ĸ��ƣ����䣩���塣
s_ref = [s_ref,zeros(Naz,Nrg-Nr)];     
% �Ը������壬��˲��㡣
 
S_ref = fft(s_ref,NFFT_rg,2);         
% ��������ľ��븵��Ҷ�任��
H_range = conj(S_ref);                 % ������ƥ���˲�����
S_range_c = S_range.*H_range;  % ����ƥ���˲�����          
% Ϊ�����Ƚϣ��ȱ任�ض�άʱ��
s_rc = ifft(S_range_c,[],2);   

% ��ͼ��ʾ
figure;
imagesc(abs(s_rc));
title('����ѹ����/��άʱ��'); 
colormap(gray);
colorbar;

%%
% �ò�����ɵĹ���Ϊ���ξ���ѹ����SRC��
% ����ѹ����Ҫ��Ƶ�����
% �ڶ�άʱ��������ݰ���
s_rc = s_rc.*exp(-1j*2*pi*f_doppler.*(ta.'*ones(1,Nrg))); 
% ��λ����Ҷ�任���������������
S_2df = fft(s_rc,NFFT_az,1);
% ��������Ҷ�任������άƵ��
S_2df = fft(S_2df,Nrg,2);               

% ��б�ӽ��µ��㶯����
D_fn_Vr = sqrt(1-lamda^2.*(fa.').^2./(4*Vr^2));         
K_src = 2*Vr^2*f0^3.*D_fn_Vr.^3./(c*R0*(fa.').^2);   
K_src_1 = 1./K_src;
% ���ξ���ѹ���˲���
H_src = exp(-1j*pi.*K_src_1*(fr.^2)); 
H_src = fftshift(H_src,2);   

S_2df_src = S_2df.*H_src;  
% ��ɶ��ξ���ѹ����SRC�����ص������������
S_rd = ifft(S_2df_src,[],2);  


% ��ͼ��ʾ
figure;
imagesc(abs(S_rd));
title('SRC��/�����������');
colormap(gray);
colorbar
%%
% �ò��ֵĹ���Ϊ��ɾ����㶯У��
% �����㶯У��Ӧ���ھ�������������
% ������һ�������Ѿ�����Ӷ�άƵ����Ƶ������������
% ÿһ�����б�ࣨR0�������ž����ŵĲ�ͬ���ı䡣
tr_RCMC = T_start + ( -Nrg/2 : (Nrg/2-1) )/Fr;   % ����RCMC�ͷ�λMF�ľ���ʱ���ᡣ
R0_RCMC = (c/2).*tr_RCMC;   % ������߱仯��R0����ΪR0_RCMC����������RCM��Ka��
delta_Rrd_fn = ((1-D_fn_Vr)./D_fn_Vr)*R0_RCMC;      % ��б�ӽ��µ�RCM

num_range = c/(2*Fr);   % һ�����������Ԫ����Ӧ�ĳ���
delta_Rrd_fn_num = delta_Rrd_fn./num_range; % ÿһ����λ��Ƶ�ʣ���RCM��Ӧ�ľ��������Ԫ��

R = 8;  % sinc��ֵ�˳���
S_rd_rcmc = zeros(NFFT_az,Nrg); % �������RCMC���ֵ
for p = 1 : NFFT_az
    for q = 1 : Nrg   % ��ʱ������ĳ����� (Nrg-Nr+1)=Nrg        
        delta_Rrd_fn_p = delta_Rrd_fn_num(p,q);
        Rrd_fn_p = q + delta_Rrd_fn_p;
        
        Rrd_fn_p = rem(Rrd_fn_p,Nrg);  % ����RCM�ĳ��Ȼᳬ��Nrg��������������һ�¡�
        
        Rrd_fn_p_int = ceil(Rrd_fn_p);        % ceil������ȡ����
        ii = ( Rrd_fn_p-(Rrd_fn_p_int-R/2):-1:Rrd_fn_p-(Rrd_fn_p_int+R/2-1)  );        
        rcmc_sinc = sinc(ii);
        rcmc_sinc = rcmc_sinc/sum(rcmc_sinc);   % ��ֵ�˵Ĺ�һ��
               
        % ����S_rdֻ��������ȡֵ���ҷ�Χ���ޡ���˲�ֵ��Ҫ��������ȡֵ����߽����⡣
        % �����Ҳ�ȡѭ����λ��˼�룬�������ȡֵ������⡣
        if (Rrd_fn_p_int-R/2) > Nrg    % ȫ����
            variable = (Rrd_fn_p_int-R/2-Nrg:1:Rrd_fn_p_int+R/2-1-Nrg);
        else
            if (Rrd_fn_p_int+R/2-1) > Nrg    % ��������
                variable_1 = (Rrd_fn_p_int-R/2:1:Nrg);
                variable_2 = (1:1:Rrd_fn_p_int+R/2-1-Nrg);
                variable = [variable_1,variable_2];
            else
                if (Rrd_fn_p_int+R/2-1) < 1    % ȫ���磨�����ܷ�����������Ҫ���ǣ�
                    variable = (Rrd_fn_p_int-R/2+Nrg:1:Rrd_fn_p_int+R/2-1+Nrg);
                else
                    if (Rrd_fn_p_int-R/2) < 1       % ��������
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
% S_rd_rcmc ����RCMC��ľ����������Ƶ�ס�

figure;
imagesc(abs(S_rd_rcmc));
title('RCMC��/�����������');
colormap(gray);
colorbar;
%%
%�ò�����ɵĹ���Ϊ��λѹ��
fa_azimuth_MF = fa;         % ��λƵ���ᣬ���ú�RCMC�����õ�Ƶ������ͬ��
Haz = exp(1j*4*pi.*(D_fn_Vr*R0_RCMC).*f0./c);   % �Ľ��ķ�λ��MF

S_rd_c = S_rd_rcmc.*Haz;            % ����ƥ���˲���
s_ac = ifft(S_rd_c,[],1);       	% ��ɷ�λѹ�����䵽ͼ���򡣽�����

%%
% ȫ�����źŴ�������Ѿ���ɡ����������г���
% ���ö�̬��ʾ��Χ������ʾͼ��
sout = (s_ac);
sout = abs(sout)/max(max(abs(sout)));
G = 20*log10(sout+eps);                         % db��ʾ
clim = [-55,0];                                           % ��̬��ʾ��Χ
figure;
imagesc(((0:Nrg-1)+first_rg_cell)/Fr*c/2+R0,((0:Naz-1)+first_rg_line)/FA*Vr,G,clim);
axis xy;
title('RADARSAT-1����/RD�㷨����')
xlabel('����(m)')
ylabel('��λ(m)')
colormap(gray);
colorbar;






