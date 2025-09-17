% main_RGPO_whitebox_strict.m
% 严格按照论文第二章 & 2.5 节复现的白盒 RGPO 仿真框架（主脚本）
clear; close all; clc;
rng(1); % 为复现性固定随机种子

% ------------------ 论文/实验参数（按论文示例） ------------------
dt = 1;                    % s (扫描周期)
K = 15;                    % 最大干扰阶段数（论文示例）
J = 20;                    % 每个候选策略重复次数（蒙特卡洛次数，论文中有 Jz 等)
% 真实目标初始状态（CV2 基准示例，单位 m / m/s）
x0 = [50e3; 50; 55e3; 350];   % (x, vx, y, vy)
% 过程噪声（论文设定 N(0,5^2)）
Qproc_xy = (5)^2 * eye(2);

% 测量噪声（论文 2.5：sigma_r = 10 m, var(theta)=0.001752)
sigma_r = 10;
sigma_theta = sqrt(0.001752);
R_polar = diag([sigma_r^2, sigma_theta^2]);

PD = 0.99;                 % 检测概率 (论文)
gate_thresh = 9.21;        % 验证门阈 (chi2 99% for 2 DOF ~9.21)
lambda_FA = 2;             % avg false alarms per scan (示例)
w1 = 0.5; w2 = 0.5 / 3000; % 评分权重示例

% ------------------ 候选策略集合（示例：也可由 PSO 输出） ------------------
N_candidates = 8;
max_dtaus_us = 50; % 每阶段最大增量（微秒）
phi_candidates = cell(N_candidates,1);
for n=1:N_candidates
    dtaus_us = rand(1,K) * max_dtaus_us; % in microseconds
    phi_candidates{n} = dtaus_us * 1e-6;  % 转换为 seconds
end

% ------------------ 循环评估每个候选策略 ------------------
scores = zeros(N_candidates,1);
details = cell(N_candidates,1);
for n = 1:N_candidates
    fprintf('Evaluating candidate %d / %d ...\n', n, N_candidates);
    rsum = 0;
    xi_list = zeros(J,1);
    MD_list = zeros(J,1);
    for j = 1:J
        % 1) 生成真实目标轨迹 Xr (4xK)
        Xr = generateTarget_CV(x0, K, dt, Qproc_xy);

        % 2) 根据候选策略 phi 生成假目标轨迹 Xf
        Xf = generateFakeFromStrategy(Xr, phi_candidates{n}, dt);

        % 3) 生成量测序列（含虚警），返回 cell 1xK，每项 Nx2 (r,theta)
        radarPos = [0,0];
        Z_all = generateMeasurements_v2(Xr, Xf, sigma_r, sigma_theta, PD, radarPos, lambda_FA, 120e3);

        % 4) 用 EKF+PDA 跟踪（严格实现）
        params = struct();
        params.dt = dt;
        params.F = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1];
        q_ax = 5^2; q_ay = 5^2;
        Q = [ (dt^4)/4*q_ax, (dt^3)/2*q_ax, 0, 0;
              (dt^3)/2*q_ax, (dt^2)*q_ax,   0, 0;
              0, 0, (dt^4)/4*q_ay, (dt^3)/2*q_ay;
              0, 0, (dt^3)/2*q_ay, (dt^2)*q_ay ];
        params.Q = Q;
        params.R_polar = R_polar;
        params.PD = PD;
        params.gate_thresh = gate_thresh;
        params.lambda_FA = lambda_FA;
        params.V_gate = [];
        params.x0 = [Xr(1,1); Xr(2,1); Xr(3,1); Xr(4,1)];
        params.P0 = diag([1e4, 1e2, 1e4, 1e2]);

        [Xest_hist, P_hist, assoc_info] = run_EKF_PDA(Z_all, params);

        % 5) 计算 xi 和 MD (按论文式(2-23),(2-24))
        [xi, MD] = computeScore_strict(Xr, Xest_hist, assoc_info, gate_thresh);
        xi_list(j) = xi; MD_list(j) = MD;

        % 6) 组合成 r (式(2-29))
        rj = scoreCombine_strict(xi, MD, w1, w2);
        rsum = rsum + rj;
    end
    scores(n) = rsum / J;
    details{n}.xi = xi_list;
    details{n}.MD = MD_list;
    fprintf('Candidate %d -> mean r = %.6f (xi mean=%.4f, MD mean=%.2f m)\n', ...
        n, scores(n), mean(xi_list), mean(MD_list));
end

[bestScore, idxBest] = max(scores);
fprintf('\\nBest candidate index %d, score=%.6f\\n', idxBest, bestScore);

% 绘示例结果（最后一次仿真）
figure; hold on; grid on;
plot(Xr(1,:), Xr(3,:),'k-','LineWidth',1.5);
plot(Xf(1,:), Xf(3,:),'b--','LineWidth',1.2);
plot(Xest_hist(1,:), Xest_hist(3,:),'r-.','LineWidth',1.3);
legend('True','Fake','Estimate'); xlabel('x (m)'); ylabel('y (m)');
title('Trajectory: true / fake / estimate (example)');
