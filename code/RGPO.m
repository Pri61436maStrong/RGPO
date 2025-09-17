% 距离波门拖引干扰建模
% 保存GIF图片
close all; clear; clc;

%% 雷达参数设置
C = 3.0e8;                  % 光速
RF = 3.75e9;                % 雷达载频 3.75GHz
Lamda = C / RF;             % 雷达工作波长
B = 5.0e6;                  % 发射信号带宽 5MHz
T = 20.0e-6;                % 发射信号时宽 20us
K = B / T;                  % 调频斜率
Fs = 20.0e6;                % 采样频率 20MHz
Ts = 1 / Fs;                % 采样间隔
PRF = 16e3;                 % 脉冲重复频率 16kHz
PRT = 1 / PRF;              % 脉冲重复间隔
M = 64;                     % 回波脉冲数 64
T_coherent = PRT * M;       % 相参积累时间（CPI）

SampleNumber = fix(Fs * PRT);   % 一个脉冲重复间隔的采样点数
number = fix(Fs * T);           % 脉冲波形的采样点数

R_max = C * PRT / 2;            % 雷达最大不模糊距离
V_max = PRF * Lamda / 4;        % 雷达最大不模糊速度
detaR = R_max / SampleNumber;   % 每点的距离间距（最小可分辨距离为7.5km，也就是说目标的间距最好为7.5的整数倍）
R_bomen = detaR;                 % 距离波门
R_min = C / (2 * B);            % 最小可探测距离
V_bomen = (C * PRF) / (2 * M * RF); % 计算雷达每个速度波门的速度
V_total = (M - 1) * V_bomen;    % 总速度范围

r = (1:SampleNumber) * detaR;   % 距离轴，这里有问题
v = (0:M-1) * V_bomen;          % 速度轴

R_gate = 50 * R_bomen;          % 距离门宽度

%% 目标参数
SNR = -5;                       % 信干比
SJR = -5;                       % 干信比（干扰/信号功率比）
R_target = 1005;                % 目标初始距离
V_target = 80;                  % 目标速度
Fd_target = 2 * V_target / Lamda; % 目标多普勒频率

target_SigPower = 1;            % 目标信号功率
jamming_SigPower = 1 / (10^(SJR/10)); % 干扰信号功率，拖引干扰功率

N = ceil(V_target / V_total); % 速度周期数

%% 干扰机参数
t_1 = 0.5;                      % 干扰停拖时间0.5s
t_2 = 3;                        % 干扰拖引时间3s
t_3 = 0.5;                      % 干扰关闭时间0.5s

N_1 = t_1 / T_coherent;         % 停拖时间对应的CPI数
N_2 = t_2 / T_coherent;         % 拖引时间对应的CPI数
N_3 = t_3 / T_coherent;         % 关闭时间对应的CPI数
N_total = N_1 + N_2 + N_3;      % 总拖引时间对应的CPI数

V_jamming = 200;                % 干扰拖引速度

%% 产生线性调频信号
% t_number = linspace(-T/2, T/2, number);
t_number = linspace(0, T, number);
Chirp = exp(1i * pi * K * t_number.^2); % LFM信号
coeff = conj(fliplr(Chirp));             % 脉压系数

%% 创建时延参数（距离）和多普勒参数（速度）
DelayNumber_target = zeros(1, M * N_total);   % 目标时延
DelayNumber_jamming = zeros(1, M * N_total);  % 干扰时延

jj = 0;
ii = 0;
for i = 1:M*N_total
    if (i > M*N_1 && i <= M*(N_2 + N_1))
        jj = jj + 1;
    else
        ii = ii + 1;
        jj = 0;
    end
    % 干扰目标的时延
    DelayNumber_jamming(1,i) = fix(Fs * 2 * (R_target + PRT * ii * V_target + jj * PRT * V_jamming) / C);
    % 真实目标的时延
    DelayNumber_target(1,i) = fix(Fs * 2 * (R_target + i * PRT * V_target) / C);
end

FreqMove_target = exp(1i * 2 * pi * Fd_target * (0:M*SampleNumber-1) / Fs);

Fd_jamming = zeros(N_total, 1);
FreqMove_jamming = zeros(N_total, M*SampleNumber);

jj = 0;
ii = 1;
for i = 1:N_total
    if (i > N_1 && i <= N_2 + N_1)
        jj = 1;
        ii = 0;
    else
        jj = 0;
    end
    Fd_jamming(i) = 2 * (jj * V_jamming + ii * V_target) / Lamda;
    FreqMove_jamming(i,:) = exp(1i * 2 * pi * Fd_jamming(i) * (0:M*SampleNumber-1) / Fs);
end

%% 主循环
Fig = figure;
set(gcf, 'position', [100 100 1400 700]);
set(gcf, 'DefaultTextFontName', 'SimHei'); % 设置默认字体
set(gcf, 'DefaultAxesFontName', 'SimHei');

filename = 'RGPO_model_1.gif';

pc = zeros(M, SampleNumber);
jj = 1;
nn = 1;

for n = 1:N_total
    if (n > N_2 + N_1)
        jj = 0; % 干扰关闭
    end
    
    for i = 1:M
        SignalTemp = zeros(1, SampleNumber);         % 初始化，目标信号模板
        SignalTemp_jamming = zeros(1, SampleNumber); % 干扰信号模板
        Signal = zeros(1, SampleNumber);             % 目标信号
        Signal_jamming = zeros(1, SampleNumber);     % 干扰信号

        % 目标回波信号(到这里)
        % idx_target = DelayNumber_target(1 + M*(n-1)) + 1;
        % SignalTemp(idx_target:idx_target+number-1) = sqrt(target_SigPower) * Chirp;
        SignalTemp(DelayNumber_target(1+M*(n-1))+1:DelayNumber_target(1+M*(n-1))+number)=sqrt(target_SigPower)*Chirp;%一个脉冲信号(未加多普勒频移)

        Signal = SignalTemp .* FreqMove_target((i-1)*SampleNumber+1 : i*SampleNumber);

        Echo = awgn(Signal, SNR, 'measured');        % 加入信噪比为10dB的噪声，加入前预估信号
        Echo_SystemNoise = Echo - Signal;

        Signal_Noise = Signal + Echo_SystemNoise;

        % 拖引干扰信号
        % idx_jamming = DelayNumber_jamming(1 + M*(n-1)) + 1;
        % SignalTemp_jamming(idx_jamming:idx_jamming+number-1) = sqrt(jamming_SigPower) * Chirp;
        SignalTemp_jamming(DelayNumber_jamming(1+M*(n-1))+1:DelayNumber_jamming(1+M*(n-1))+number)=sqrt(jamming_SigPower)*Chirp;
        Signal_jamming = SignalTemp_jamming .* FreqMove_jamming(n, (i-1)*SampleNumber+1 : i*SampleNumber);

        Signal_jamming_Noise = Signal_jamming + Echo_SystemNoise;

        % 雷达总回波（目标 + 干扰 + 噪声）
        % jj ——当干扰机关闭时，雷达回波中没有干扰信号
        Echo = jj * Signal_jamming + Signal + Echo_SystemNoise;

        % 频域脉压
        Echo_fft = fft(Echo, SampleNumber+number-1);
        coeff_fft = fft(coeff, SampleNumber+number-1);
        pc_fft = Echo_fft .* coeff_fft;
        pc_freq0 = ifft(pc_fft);
        pc_freq1 = pc_freq0(number : SampleNumber+number-1);
        pc(i, :) = pc_freq1;
    end

    %% 绘制脉压前信号实部，雷达接收回波的时域变化图
    subplot(2,1,1);
    % plot(real(Signal_Noise), 'black', 'LineWidth', 1.5);
    plot(real(Signal), 'black', 'LineWidth', 2);
    hold on;
    % plot(real(jj * Signal_jamming_Noise), 'r', 'LineWidth', 1.5);
    plot(real(jj * Signal_jamming), 'r', 'LineWidth', 2);
    hold off;
    set(gca, 'FontSize', 18, 'FontName', 'SimHei');
    xlabel('采样点'); ylabel('幅度(V)');
    axis([0, SampleNumber, -2.5, 2.5]);

    str_1 = [num2str(n*M*PRT, '%.3f') 's（RGPO捕获期）时脉压前信号的实部'];
    str_2 = [num2str(n*M*PRT, '%.3f') 's（RGPO拖引期）时脉压前信号的实部'];
    str_3 = [num2str(n*M*PRT, '%.3f') 's（RGPO停止期）时脉压前信号的实部'];
    
    if (n <= N_1)
        title(cellstr(str_1));
    elseif n > N_1 && n <= N_2+N_1
        title(cellstr(str_2));
    else
        title(cellstr(str_3));
    end
    box on;
    legend('目标回波', '干扰回波');
    grid on;

    %% 绘制脉压后结果（距离维实时画图）
    Z = abs(pc(M/2, :));
    Z = Z / max(Z);
    Z = 20 * log10(Z);
    if (n > N_2 + N_1) % 干扰关闭时，只有目标回波
        Z = Z + SJR;
    end
    [PKS, LOCS] = findpeaks(Z, 'NPeaks', 1, 'SortStr', 'descend');
    
    % if isempty(LOCS)
    %     LOCS = 1;
    % end

    if (n > N_1 + N_2 + 1)
        nn = nn + 1;
        L(nn + N_1 + N_2) = LOCS;
        if L(nn + N_1 + N_2) < L(nn + N_1 + N_2 - 1)
            L(nn + N_1 + N_2) = L(nn + N_1 + N_2 - 1);
        end
    elseif (n <= N_1 + N_2 + 1)
        L(n) = LOCS;
    end

    str = [num2str((L(n)-1)*detaR) 'm'];

    % 设置距离波门
    if (n <= N_1 + N_2)
        LOCS_ref = L(n);
    end
    Gate = -60 * ones(1, SampleNumber);
    Gate(LOCS_ref-(R_gate/2)/R_bomen:LOCS_ref+(R_gate/2)/R_bomen)=0;
    % gate_indices = round(LOCS_ref - (R_gate/2)/detaR : LOCS_ref + (R_gate/2)/detaR);
    % gate_indices = gate_indices(gate_indices > 0 & gate_indices <= SampleNumber);
    % Gate(gate_indices) = 0;

    subplot(2,1,2);
    plot(r, Z, 'black', 'LineWidth', 2);
    hold on;
    plot(r, Gate, 'red', 'LineWidth', 3);
    hold off;
    text((L(n)-1)*detaR, PKS+3, cellstr(str), 'FontSize', 13);
    axis([0, 3000, -45, 6]);
    xlabel('距离（m）'); ylabel('幅度（dB）');
    set(gca, 'FontSize', 18, 'FontName', 'SimHei');
    grid on;

    str_1 = [num2str(n*M*PRT, '%.3f') 's（RGPO捕获期）时脉压后信号'];
    str_2 = [num2str(n*M*PRT, '%.3f') 's（RGPO拖引期）时脉压后信号'];
    str_3 = [num2str(n*M*PRT, '%.3f') 's（RGPO停止期）时脉压后信号'];
    
    if n <= N_1
        title(cellstr(str_1));
    elseif n > N_1 && n <= N_2+N_1
        title(cellstr(str_2));
    else
        title(cellstr(str_3));
    end
    box on;
    legend('脉压处理结果', '距离波门');
    drawnow;

    %% 保存GIF
    frame = getframe(Fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 1024);
    if n == 1
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'overwrite', 'Loopcount', inf);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
    disp(n);
end