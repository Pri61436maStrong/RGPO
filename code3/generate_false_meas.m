%% generate_false_meas.m：生成虚假目标状态及量测（含量测噪声）
% 输入：Xr - 真实目标状态序列（维度×K），phi - 干扰策略（Δτ序列，1×K）
% 输出：Zf - 虚假目标量测序列（2×K，[距离; 角度]）
function [Xf,Zf] = generate_false_meas(Xr, phi)
    global RadarParams ExpParams;
    c = RadarParams.c;
    % T = RadarParams.T;
    K = ExpParams.K;
    sigma_r = RadarParams.sigma_r;         % 距离噪声标准差（10m，表2-1）
    sigma_theta = RadarParams.sigma_theta; % 角度噪声标准差（0.00175rad，表2-1）

    dim_r = size(Xr,1);  % 真实目标状态维度
    
    % 1. 初始化虚假目标状态（与真实目标同维度，仅位置有偏移）
    Xf = zeros(dim_r, K);
    Xf(:,1) = Xr(:,1);  % 初始状态与真实目标一致
    theta_k = zeros(1, K);      % 修正：初始化1×K序列，存储各阶段方位角

    % 2. 计算各阶段转发时延τk（Δτ为策略输入，τk=sum(Δτ_1~Δτ_k)）
    tau = cumsum(phi);  % 转发时延序列累积和（s），其中第i个元素的值等于phi中第 1 到第i个元素的总和（即累积和）。
    Zf =  zeros(2, K); % 初始化虚假目标量测（含噪声）
    for k = 1:K
        % 2.1 真实目标位置（x_r, y_r）
        x_r = Xr(1,k);
        y_r = Xr(3,k);
        theta_k(k) = atan2(y_r, x_r);  % 逐阶段存储方位角，符合论文动态偏移设定,真实目标相对于雷达的方位角

        % 2.2 虚假目标与真实目标距离dk = c*tau(k)/2（论文1.17）
        dk = c * tau(k) / 2;

        % % 运动方向基于方位角的坐标变换，对应公式2-6
        x_f = x_r + dk * cos(theta_k(k)); %虚假目标x坐标（状态第1行）
        y_f = y_r + dk * sin(theta_k(k));
        
        % 2.4 虚假目标状态（公式2-6）
        Xf(1,k) = x_f; 
        Xf(3,k) = y_f;
        if dim_r >=4  % 速度维度噪声
            Xf(2,k) = Xr(2,k) ;  % 速度保持不变
            Xf(4,k) = Xr(4,k); 
        end
        if dim_r ==6  % CA的加速度维度
            Xf(5,k) = Xr(5,k) ;
            Xf(6,k) = Xr(6,k) ;
        end
        if dim_r ==5  % CT的角速度维度
            Xf(5,k) = Xr(5,k) ;
        end
        
        % 3. 量测模型（距离+角度，加量测噪声，论文1.25-1.33）
        
        % 2. 计算理想量测（无噪声的距离、角度）
        f_ideal = sqrt(x_f^2 + y_f^2);          % 理想距离（雷达到虚假目标的直线距离）
        f_theta_ideal = atan2(y_f, x_f);          % 理想角度（虚假目标相对于雷达的方位角）
        
        % 3. 叠加表2-1的量测值，距离和角度噪声（高斯分布）
        f_noisy =  f_ideal + normrnd(0, sigma_r);        % 距离量测+噪声（N(0,10²)）
        f_theta_noisy = f_theta_ideal + normrnd(0, sigma_theta);  % 角度量测+噪声（N(0,0.00175²)）
        
        % 6. 极坐标转换为x/y轴量测误差（即公式2-10的Wf(k)，文档1-32：δf,x/δf,y）
        Zf(1,k) = f_noisy *  cos(f_theta_noisy);  % 第一行：x
        Zf(2,k) = f_noisy * sin(f_theta_noisy);  %  第二行：y
    end
end