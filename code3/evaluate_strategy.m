%% 子函数：评估单个干扰策略的性能（J次模拟）
function [track_break_rate, avg_md, score] = evaluate_strategy(phi, benchmark_type, J)
    global ExpParams;
    K = ExpParams.K;
    w1 = ExpParams.perf_w1;
    w2 = ExpParams.perf_w2;
    track_break = 0;
    md_sum = 0;

    for j =1:J
        % 1. 生成真实目标状态及量测
        [Xr,Zr] = generate_true_target(benchmark_type); % 对应公式2-3

        % 3. 生成虚假目标状态及量测
        [Xf,Zf] = generate_false_meas(Xr, phi);

        % 4. 雷达跟踪
        [X_est, Z_valid] = radar_tracking(Zr, Zf, benchmark_type);

        % 5. 计算当前模拟的航迹中断与误差距离（论文1.75-1.85）
        z_valid_K = Z_valid{K};

        is_break = 1;
        if ~isempty(z_valid_K)
            for i =1:size(z_valid_K,2)
                if norm(z_valid_K(:,i) - Zr(:,K)) < 1e-6  % 量测匹配（误差<1e-3），z_valid_K只要一个满足，就代表没丢失目标
                    is_break = 0;
                    break;
                end
            end
        end
        track_break = track_break + is_break;

        % 5.2 误差距离：估计位置与真实位置的欧氏距离（论文1.84）
        x_est = X_est(1,K);
        y_est = X_est(3,K);
        x_r = Xr(1,K);
        y_r = Xr(3,K);
        md = sqrt((x_est - x_r)^2 + (y_est - y_r)^2);
        md_sum = md_sum + md;
    end

    % 6. 性能指标计算（正确，无修改）
    track_break_rate = track_break / J;  % 航迹中断率（期望）
    avg_md = md_sum / J;                 % 平均误差距离（期望）
    score = w1 * track_break_rate + w2 * avg_md;  % 综合得分
end
