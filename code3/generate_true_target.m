function [Xr, Zr] = generate_true_target(benchmark_type)
% 返回Xr真实目标状态矩阵
    global RadarParams PSOParams ExpParams Benchmark;

    T = RadarParams.T;
    K = ExpParams.K;
    X0 = Benchmark.(benchmark_type).X0;  % 目标初始状态[x(km)→m, ẋ(m/s), y(km)→m, ẏ(m/s)]
    Qp = Benchmark.(benchmark_type).Qp;
    dim = length(X0);  % 状态维度（CV:4, CA:6, CT:5, ERV:4）
    Xr = zeros(dim, K);
    Xr(:,1) = X0;      % 初始状态
    R = Benchmark.(benchmark_type).R;

    for k = 2:K
        switch benchmark_type
            case 'CV'  % 匀速运动（修复后）
                F = [1, T, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, T;
                     0, 0, 0, 1];
                G = Benchmark.(benchmark_type).G ;
                Xr(:,k) = F * Xr(:,k-1) + G*mvnrnd(zeros(2,1), R)';% 包含过程噪声，公式2-4的U（t-1）                
        end

    end
        % 量测矩阵 H
    H = RadarParams.H;

%% -------------------------------
    sigma_r = RadarParams.sigma_r;         % 距离噪声标准差（10m，表2-1）
    sigma_theta = RadarParams.sigma_theta; % 角度噪声标准差（0.00175rad，表2-1）
    % 新增：生成每个阶段的真实目标量测（含噪声）
    Zr0 = H * Xr; % 初始化真实目标量测
    Zr =  zeros(2, K); % 初始化真实目标量测（含噪声）
    for k = 1:K
        % 1. 提取真实目标当前阶段的位置（x_r, y_r）
        x_r = Zr0(1, k);  % 真实目标x坐标
        y_r = Zr0(2, k);  % 真实目标y坐标

        % 2. 计算理想量测（无噪声的距离、角度）
        r_ideal = sqrt(x_r^2 + y_r^2);          % 理想距离（雷达到真实目标的直线距离）
        theta_ideal = atan2(y_r, x_r);          % 理想角度（真实目标相对于雷达的方位角）

        % 3. 叠加表2-1的量测噪声（高斯分布）
        r_noisy =  r_ideal + normrnd(0, sigma_r);        % 距离量测+噪声（N(0,10²)）
        theta_noisy = theta_ideal + normrnd(0, sigma_theta);  % 角度量测+噪声（N(0,0.00175²)）

        % 6. 极坐标转换为x/y轴量测误差（即公式2-9的Wf(k)，文档1-32：δf,x/δf,y）
        Zr(1,k) = r_noisy *  cos(theta_noisy);  % 第一行：x
        Zr(2,k) = r_noisy * sin(theta_noisy);  %  第二行：y
    end
end