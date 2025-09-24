%% compare_rgpo.m：对比算法（UV-RGPO/UA-RGPO）性能评估
% 输入：benchmark_type - 基准问题类型
% 输出：uv_best - UV-RGPO最优策略及性能，ua_best - UA-RGPO最优策略及性能
function [uv_best, ua_best] = compare_rgpo(benchmark_type)
    global CompParams ExpParams PSOParams RadarParams;
    vpo_range = CompParams.UV.vpo_range;
    vipo_range = CompParams.UA.vipo_range;
    apo_range = CompParams.UA.apo_range;
    K = ExpParams.K;
    J = PSOParams.J;  % 同PSO的模拟预算（公平对比）
    w1 = ExpParams.perf_w1;
    w2 = ExpParams.perf_w2;
    
    %% 1. UV-RGPO（匀速拖引，论文1.163）
    uv_perf = zeros(length(vpo_range), 4);  % [vpo, 航迹中断率, 平均误差距离, 综合得分]
    for i =1:length(vpo_range)
        vpo = vpo_range(i);
        % 生成UV-RGPO策略：Δτ_k = 2*vpo*T/c（匀速拖引，距离增量=vpo*T）
        delta_tau = 2 * vpo * ExpParams.T / RadarParams.c * ones(1,K);
        % 评估性能
        [track_break, avg_md, score] = evaluate_strategy(delta_tau, benchmark_type, J);
        uv_perf(i,:) = [vpo, track_break, avg_md, score];
    end
    % 最优UV-RGPO（综合得分最大）
    [~, uv_idx] = max(uv_perf(:,4));
    uv_best.vpo = uv_perf(uv_idx,1);
    uv_best.phi = 2 * uv_best.vpo * ExpParams.T / RadarParams.c * ones(1,K);
    uv_best.track_break = uv_perf(uv_idx,2)/J;
    uv_best.avg_md = uv_perf(uv_idx,3)/J;
    uv_best.score = uv_perf(uv_idx,4);
    
    %% 2. UA-RGPO（匀加速拖引，论文1.163）
    ua_perf = zeros(length(vipo_range)*length(apo_range), 5);  % [vipo, apo, ...]
    cnt =1;
    for i =1:length(vipo_range)
        vipo = vipo_range(i);
        for j =1:length(apo_range)
            apo = apo_range(j);
            % 生成UA-RGPO策略：Δτ_k = 2*(vipo*T + 0.5*apo*T²)/c（加速拖引）
            delta_tau = 2 * (vipo * ExpParams.T + 0.5 * apo * ExpParams.T^2) / RadarParams.c * ones(1,K);
            % 评估性能
            [track_break, avg_md, score] = evaluate_strategy(delta_tau, benchmark_type, J);
            ua_perf(cnt,:) = [vipo, apo, track_break, avg_md, score];
            cnt = cnt +1;
        end
    end
    % 最优UA-RGPO（综合得分最大）
    [~, ua_idx] = max(ua_perf(:,5));
    ua_best.vipo = ua_perf(ua_idx,1);
    ua_best.apo = ua_perf(ua_idx,2);
    ua_best.phi = 2 * (ua_best.vipo * ExpParams.T + 0.5 * ua_best.apo * ExpParams.T^2) / RadarParams.c * ones(1,K);
    ua_best.track_break = ua_perf(ua_idx,3)/J;
    ua_best.avg_md = ua_perf(ua_idx,4)/J;
    ua_best.score = ua_perf(ua_idx,5);
    
    fprintf('UV-RGPO最优：vpo=%.1f m/s，综合得分=%.4f\n', uv_best.vpo, uv_best.score);
    fprintf('UA-RGPO最优：vipo=%.1f m/s, apo=%.1f m/s²，综合得分=%.4f\n', ua_best.vipo, ua_best.apo, ua_best.score);
end