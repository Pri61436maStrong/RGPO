%% main_exp.m：实验主函数
clear; clc; close all;
run('params.m');  % 加载参数

benchmark_list = {'CV'};
for b_idx =1:length(benchmark_list)
    benchmark_type = benchmark_list{b_idx};
    fprintf('\n==================== %s基准问题实验 ====================\n', benchmark_type);
    
    %% 2. 运行PSO-ER白盒RGPO
    fprintf('\n1. 运行PSO-ER白盒RGPO优化...\n');
    [pso_phi, pso_perf] = pso_er_rgpo(benchmark_type);
    
    %% 3. 运行对比算法（UV-RGPO/UA-RGPO）
    fprintf('\n2. 运行对比算法（UV-RGPO/UA-RGPO）...\n');
    [uv_best, ua_best] = compare_rgpo(benchmark_type);
    
    %% 4. 多次独立重复实验（降低随机误差，论文1.159）
    fprintf('\n3. %d次独立重复实验，计算平均性能...\n', ExpParams.repeat);
    perf_summary = zeros(3, 3);  % 行：PSO/UV/UA；列：航迹中断率/平均误差距离/综合得分
    for r =1:ExpParams.repeat
        % PSO-ER性能
        [pso_tb, pso_md, pso_sc] = evaluate_strategy(pso_phi, benchmark_type, PSOParams.J);
        % UV-RGPO性能
        [uv_tb, uv_md, uv_sc] = evaluate_strategy(uv_best.phi, benchmark_type, PSOParams.J);
        % UA-RGPO性能
        [ua_tb, ua_md, ua_sc] = evaluate_strategy(ua_best.phi, benchmark_type, PSOParams.J);
        
        perf_summary(1,:) = perf_summary(1,:) + [pso_tb, pso_md, pso_sc];
        perf_summary(2,:) = perf_summary(2,:) + [uv_tb, uv_md, uv_sc];
        perf_summary(3,:) = perf_summary(3,:) + [ua_tb, ua_md, ua_sc];
    end
    perf_summary = perf_summary / ExpParams.repeat;  % 平均性能

    %% 5. 输出性能表格（对应论文表2-2/2-3/2-4）
    fprintf('\n==================== %s基准问题性能对比 ====================\n', benchmark_type);
    % 首先检查perf_summary的维度
fprintf('perf_summary的维度: %d行, %d列\n', size(perf_summary));
    perf_table = table(...
        {'PSO-ER RGPO'; 'UV-RGPO'; 'UA-RGPO'}, ...
        perf_summary(:,1)*100, ...
        perf_summary(:,2), ...
        perf_summary(:,3), ...
        'VariableNames', {'策略', '航迹中断率(%)', '平均误差距离(m)', '综合干扰性能'});%
    disp(perf_table);
    
    %% 6. 绘制性能对比图（对应论文图2-7/2-8/2-9/2-10）
    figure('Position', [100, 100, 1000, 400]);
    % 子图1：航迹中断率对比
    subplot(1,3,1);
    bar([1,2,3], perf_summary(:,1)*100, 'FaceColor', [0.2,0.6,0.8]);
    set(gca, 'XTickLabel', {'PSO-ER RGPO', 'UV-RGPO', 'UA-RGPO'});
    ylabel('航迹中断率(%)');
    title([benchmark_type, '基准问题：航迹中断率对比']);
    grid on;
    
    % 子图2：平均误差距离对比
    subplot(1,3,2);
    bar([1,2,3], perf_summary(:,2), 'FaceColor', [0.8,0.4,0.2]);
    set(gca, 'XTickLabel', {'PSO-ER RGPO', 'UV-RGPO', 'UA-RGPO'});
    ylabel('平均误差距离(m)');
    title([benchmark_type, '基准问题：平均误差距离对比']);
    grid on;
    
    % 子图3：综合性能对比
    subplot(1,3,3);
    bar([1,2,3], perf_summary(:,3), 'FaceColor', [0.4,0.8,0.2]);
    set(gca, 'XTickLabel', {'PSO-ER RGPO', 'UV-RGPO', 'UA-RGPO'});
    ylabel('综合干扰性能');
    title([benchmark_type, '基准问题：综合性能对比']);
    grid on;
end