function [Z_all] = generateMeasurements_v2(Xr, Xf, sigma_r, sigma_theta, PD, radarPos, lambda_FA, r_max)
% 生成包含真目标、假目标和虚警量测（每时刻 Nx2: [r,theta]）
K = size(Xr,2);
Z_all = cell(1,K);
for k = 1:K
    Zk = [];
    if rand <= PD
        xr = Xr(1,k) - radarPos(1);
        yr = Xr(3,k) - radarPos(2);
        [r_true, theta_true] = cart2pol(xr, yr);
        r_meas = r_true + sigma_r * randn;
        theta_meas = theta_true + sigma_theta * randn;
        Zk = [Zk; r_meas, theta_meas];
    end
    if rand <= PD
        xf = Xf(1,k) - radarPos(1);
        yf = Xf(3,k) - radarPos(2);
        [r_fake, theta_fake] = cart2pol(xf, yf);
        r_meas_f = r_fake + sigma_r * randn;
        theta_meas_f = theta_fake + sigma_theta * randn;
        Zk = [Zk; r_meas_f, theta_meas_f];
    end
    N_FA = poissrnd(lambda_FA);
    if N_FA > 0
        r_FA = r_max * rand(N_FA,1);
        theta_FA = -pi + 2*pi*rand(N_FA,1);
        Zk = [Zk; [r_FA, theta_FA]];
    end
    Z_all{k} = Zk;
end
end
