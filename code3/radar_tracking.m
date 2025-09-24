%% radar_tracking.m：雷达跟踪（KF-PDA/EKF-PDA），输出跟踪结果
% 输入：Zr - 真实目标量测，Zf - 虚假目标量测，benchmark_type - 基准问题类型
% 输出：X_est - 目标状态估计序列（维度×K），Z_valid - 各阶段有效量测集
function [X_est, Z_valid] = radar_tracking(Zr, Zf, benchmark_type)
    global RadarParams ExpParams Benchmark;
    T = RadarParams.T;
    K = ExpParams.K;
    gamma = RadarParams.gamma; % 跟踪波门阈值
    P_D = RadarParams.P_D; % 检测概率

    tracking_type = Benchmark.(benchmark_type).tracking;
    X0 = Benchmark.(benchmark_type).X0;
    Qp = Benchmark.(benchmark_type).Qp;
    Qm = RadarParams.Qm;
    dim = length(X0);
    H = RadarParams.H;
    
    % 1. 初始化跟踪状态
    X_est = zeros(dim, K);
    X_est(:,1) = X0;  % 初始估计=真实初始状态
    P_est = eye(dim) * 1e4;  % 初始协方差（大值，代表不确定性）
    Z_valid = cell(1, K);    % 有效量测集（每个阶段存多个量测）   

    for k = 2:K
        % 2. 状态预测（论文1.54）
        if strcmp(tracking_type, 'KF-PDA')
            % KF：线性预测
            switch benchmark_type
                case 'CV'
                    F = [1, T, 0, 0; 0, 1, 0, 0; 0, 0, 1, T; 0, 0, 0, 1]; % 状态转移矩阵，对应2-16
            end
            X_pred = F * X_est(:,k-1);  % 一步预测，公式2-14
            Z_pred = H * X_pred;  % 量测预测，对应2-15
            P_pred = F * P_est * F' + Qp;  % 协方差预测,公式2-16，Q为表格2-1中的100        
        end     
              
        S = H * P_pred * H' + Qm;  % 新息协方差    2-18
        S_inv = inv(S);
        % 4. 有效量测筛选（波门内的量测：(Z-Z_pred)^T*S_inv*(Z-Z_pred) ≤ gamma）
        Z_all = [Zr(:,k), Zf(:,k)];  % 真实+虚假量测,Zc候选量测，2-17
        Z_valid_k = []; % 有效量测

        for i = 1:size(Z_all,2)
            delta = Z_all(:,i) - Z_pred;
            if delta' * S_inv * delta <= gamma  % 公式2-17
                Z_valid_k = [Z_valid_k, Z_all(:,i)];
            end
        end
        Z_valid{k} = Z_valid_k;
        zeta = size(Z_valid_k,2);  % 有效量测个数
        
        % 5. PDA关联概率计算（论文1.58）
        if zeta ==0
            rho = 1;  % 无有效量测，关联概率=1（无目标）
            X_update = X_pred;
            P_update = P_pred;
        else
            % 每个有效量测的似然度
            lik = zeros(1, zeta);
            for i =1:zeta
                delta = Z_valid_k(:,i) - Z_pred;
                lik(i) = (1/(2*pi*sqrt(det(S)))) * exp(-0.5*delta'*S_inv*delta);
            end
            % 关联概率（含“无目标”项rho0）
            rho0 = (1 - P_D) / (P_D * zeta + (1 - P_D));
            rho = (P_D * lik) / (P_D * zeta + (1 - P_D));
            rho = [rho0, rho];  % rho(1)=rho0，rho(2:zeta+1)=各量测关联概率
        end
        
        % 6. 状态更新（论文1.63）
        X_update = X_pred;
        P_update = P_pred;
        if zeta >0
            for i =1:zeta
                delta = Z_valid_k(:,i) - Z_pred;
                if strcmp(tracking_type, 'KF-PDA')
                    K_gain = P_pred * H' * S_inv;  % KF卡尔曼增益
                else
                    K_gain = P_pred * H' * S_inv;  % EKF增益（同KF形式）
                end
                X_i = X_pred + K_gain * delta;  % 第i个量测的状态更新
                P_i = (eye(dim) - K_gain * H) * P_pred;  % 第i个量测的协方差更新
                % 加权更新（含关联概率）
                X_update = X_update + rho(i+1) * (X_i - X_pred);
                P_update = P_update + rho(i+1) * (P_i - P_pred + (X_i - X_update)*(X_i - X_update)');
            end
        end
        X_est(:,k) = X_update;
        P_est = P_update;

    end
end