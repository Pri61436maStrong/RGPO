%% params.m：实验核心参数（对应论文表2-1、1.159）
clear; clc; close all;

%% 1. 雷达通用参数（优先定义，因后续基准问题依赖此参数）
global RadarParams;
RadarParams.T = 1;                % 扫描间隔时间（s），论文1.154（表2-1）
RadarParams.sigma_r = 10;         % 距离量测噪声标准差（m），N(0,10²)，论文1.154
RadarParams.sigma_theta = 0.00175;% 角度量测噪声标准差（rad），N(0,0.00175²)，论文1.154
RadarParams.P_D = 0.99;           % 检测概率，论文1.154
RadarParams.gamma = 16;           % 跟踪波门阈值，论文1.154
RadarParams.c = 3e8;              % 光速（m/s），论文1.17
RadarParams.H = [1 0 0 0; 0 0 1 0]; % 量测矩阵，跟踪雷达，只保留位置信息
RadarParams.R = diag([RadarParams.sigma_r^2, RadarParams.sigma_theta^2]);
RadarParams.Qm = diag([RadarParams.sigma_r^2, RadarParams.sigma_theta^2]);

%% 2. 算法参数（PSO-ER）
global PSOParams;
PSOParams.N = 10;                  % 粒子数40，论文1.159
PSOParams.J = 2000;                % 每个候选策略的模拟预算2000，论文1.159
PSOParams.gmax = 100;              % 最大迭代次数100，论文1.159
PSOParams.c1 = 1.49445;            % 个体学习因子，论文1.137
PSOParams.c2 = 1.49445;            % 全局学习因子，论文1.137
PSOParams.w_min = 0.4;             % 最小惯性权重，论文1.137
PSOParams.w_max = 0.9;             % 最大惯性权重，论文1.137

%% 3. 基准问题参数（CV/CA/CT/ERV，对应论文1.155-1.158）
% 关键修复：所有依赖“扫描间隔时间T”的计算，均引用 RadarParams.T
global Benchmark;
% 3.1 CV基准问题（匀速运动）
Benchmark.CV.X0 = [50e3; 50; 55e3; 350];  % 初始状态[x(km)→m, ẋ(m/s), y(km)→m, ẏ(m/s)]
Benchmark.CV.G = [0.5*RadarParams.T^2, 0; 
                 RadarParams.T,       0; 
                 0, 0.5*RadarParams.T^2; 
                 0,       RadarParams.T];
Benchmark.CV.R = diag([5^2,5^2]);  % 目标过程噪声，G*R对应公式2-4的U（t-1）
Benchmark.CV.v_min = 2;         % m/s，论文2.29
Benchmark.CV.v_max = 50;           % m/s
Benchmark.CV.Qp = [25,50,0,0;50,100,0,0;0,0,25,50;0,0,50,100];           % 雷达跟踪方法为 KF-PDA,x、y 轴上的过程噪声协方差均设置为 100。
% Benchmark.CV.Qp = diag([100,100,100,100]);        % 雷达跟踪方法为 KF-PDA,x、y 轴上的过程噪声协方差均设置为 100。
Benchmark.CV.tracking = 'KF-PDA';         % 跟踪算法，论文1.155

%% 4. 对比算法参数（UV-RGPO/UA-RGPO，论文1.163）
global CompParams;
CompParams.UV.vpo_range = 2:2:50;    % UV-RGPO拖引速度遍历范围（m/s）
CompParams.UA.vipo_range = 2:2:50;   % UA-RGPO拖引初速度范围
CompParams.UA.apo_range = 2:2:20;    % UA-RGPO拖引加速度范围

%% 5. 实验控制参数
global ExpParams;
ExpParams.K = 15;                     % 最大干扰阶段数（可设10/15/20，论文1.164）
ExpParams.repeat = 1;                % 独立重复实验次数50，论文1.159
ExpParams.perf_w1 = 0.5;              % 航迹中断率权重，论文1.164
ExpParams.perf_w2 = 0.5/200;          % 平均误差距离权重，论文1.164
ExpParams.T = 1;