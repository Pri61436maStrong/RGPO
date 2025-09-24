%% pso_er_rgpo.m：PSO-ER优化的白盒RGPO干扰策略
% 输入：benchmark_type - 基准问题类型
% 输出：phi_opt - 最优干扰策略（1×K），perf_opt - 最优性能（航迹中断率,平均误差距离,综合得分）
function [phi_opt, perf_opt] = pso_er_rgpo(benchmark_type)
    global PSOParams ExpParams Benchmark RadarParams;
    N = PSOParams.N;    % 粒子数
    J = PSOParams.J;    % 模拟预算
    gmax = PSOParams.gmax;  % 最大迭代
    K = ExpParams.K;    % 干扰阶段数
    w1 = ExpParams.perf_w1;
    w2 = ExpParams.perf_w2;


    % 1. 初始化粒子（拉丁超立方抽样LHS）
    % Δτ范围：13ns~333ns（转化为秒）

    tau_min =  2*Benchmark.(benchmark_type).v_min*RadarParams.T/RadarParams.c;
    tau_max = 2*Benchmark.(benchmark_type).v_max*RadarParams.T/RadarParams.c;
    lhs_sample = lhsdesign(N, K);  % LHS抽样（N行K列，[0,1]范围）
    phi = tau_min + lhs_sample * (tau_max - tau_min);  % 映射到Δτ范围（粒子位置）,phi就是时延策略，每一行表示一个粒子的Δτ
    v = 0.1 * (tau_max - tau_min) * (rand(N,K) - 0.5);  % 粒子速度（初始小范围）
    
    % 2. 初始化个体最优与全局最优
    pbest_phi = phi;  % 个体最优位置（初始=粒子位置）
    pbest_perf = zeros(N, 3);  % 个体最优性能（航迹中断率,平均误差距离,综合得分）
    for n =1:N
        % 每个粒子的初始性能评估（J次模拟）
        [track_break, avg_md, score] = evaluate_strategy(phi(n,:), benchmark_type, J);
        pbest_perf(n,:) = [track_break, avg_md, score];
    end
    [~, g_idx] = max(pbest_perf(:,3));  % 全局最优索引（综合得分最大）
    gbest_phi = pbest_phi(g_idx,:);
    gbest_perf = pbest_perf(g_idx,:);
    
    % 3. PSO迭代优化
    for g =1:gmax
        % 3.1 计算惯性权重（线性递减：0.9→0.4，论文1.137）
        w = PSOParams.w_max - (PSOParams.w_max - PSOParams.w_min) * g/gmax; % 对应公式2-35
        
        % 3.2 更新粒子速度与位置（论文1.137-1.138）
        delta1 = rand(N,K);  % [0,1]随机数
        delta2 = rand(N,K);
        v = w * v + PSOParams.c1 * delta1 .* (pbest_phi - phi) + PSOParams.c2 * delta2 .* (gbest_phi - phi);
        % 速度限制（避免溢出，±20%Δτ范围）
        v = min(v, 0.2*(tau_max - tau_min));
        v = max(v, -0.2*(tau_max - tau_min));
        phi = phi + v;
        % 位置限制（Δτ不能超出[tau_min, tau_max]）
        phi = min(phi, tau_max);
        phi = max(phi, tau_min);
        
        % 3.3 评估当前粒子性能（ER：每个粒子J次模拟，论文1.131）
        curr_perf = zeros(N,3);
        for n =1:N
            [track_break, avg_md, score] = evaluate_strategy(phi(n,:), benchmark_type, J);
            curr_perf(n,:) = [track_break, avg_md, score];
        end
        
        % 3.4 更新个体最优与全局最优
        for n =1:N
            if curr_perf(n,3) > pbest_perf(n,3)  % 综合得分更高则更新
                pbest_phi(n,:) = phi(n,:);
                pbest_perf(n,:) = curr_perf(n,:);
            end
        end
        [~, g_idx_new] = max(pbest_perf(:,3));
        if pbest_perf(g_idx_new,3) > gbest_perf(3)
            gbest_phi = pbest_phi(g_idx_new,:);
            gbest_perf = pbest_perf(g_idx_new,:);
        end
        
        % 3.5 迭代日志
        fprintf('PSO迭代%d/%d，当前全局最优综合得分：%.4f\n', g, gmax, gbest_perf(3));
    end
    
    % 4. 输出最优策略与性能
    phi_opt = gbest_phi;
    perf_opt = gbest_perf;
end

