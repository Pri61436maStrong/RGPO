function [Xest_hist, P_hist, assoc_info] = run_EKF_PDA(Z_all, params)
% EKF measurement model (polar) + PDA association (practical implementation)

dt = params.dt;
F = params.F;
Q = params.Q;
R_polar = params.R_polar;
PD = params.PD;
gate_thresh = params.gate_thresh;
lambda_c = [];
use_lambda_FA = false;
if isfield(params,'lambda_c') && ~isempty(params.lambda_c)
    lambda_c = params.lambda_c;
elseif isfield(params,'lambda_FA') && isfield(params,'V_gate') && ~isempty(params.lambda_FA) && ~isempty(params.V_gate)
    lambda_FA = params.lambda_FA;
    V_gate = params.V_gate;
    use_lambda_FA = true;
end

K = numel(Z_all);
if isfield(params,'x0'), x = params.x0; else x = zeros(4,1); end
if isfield(params,'P0'), P = params.P0; else P = diag([1e6,1e6,1e6,1e6]); end

Xest_hist = zeros(4,K);
P_hist = zeros(4,4,K);
assoc_info = struct();
assoc_info.beta = cell(1,K);
assoc_info.numMeas = zeros(1,K);
assoc_info.likelihoods = cell(1,K);

I4 = eye(4);

for k = 1:K
    % Predict
    xpred = F * x;
    Ppred = F * P * F' + Q;

    Zk = Z_all{k}; % Nx2 [r theta]
    if isempty(Zk)
        x = xpred;
        P = Ppred;
        Xest_hist(:,k) = x;
        P_hist(:,:,k) = P;
        assoc_info.beta{k} = [];
        assoc_info.numMeas(k) = 0;
        assoc_info.likelihoods{k} = [];
        continue;
    end

    % For each measurement compute innovation, H, S, likelihood
    Nmeas = size(Zk,1);
    nu = zeros(2,Nmeas);
    Hi_cell = cell(1,Nmeas);
    Si = zeros(2,2,Nmeas);
    L = zeros(Nmeas,1);

    px = xpred(1); py = xpred(3);
    rpred = sqrt(px^2 + py^2);
    thetapred = atan2(py, px);
    zpred = [rpred; thetapred];

    for i = 1:Nmeas
        z = Zk(i,:)';
        dtheta = wrapToPi(z(2) - zpred(2));
        nu(:,i) = [ z(1) - zpred(1); dtheta ];

        if rpred < 1e-6
            rpred = 1e-6;
        end
        Hi = zeros(2,4);
        Hi(1,1) = px / rpred;  Hi(1,3) = py / rpred;
        Hi(2,1) = -py / (rpred^2); Hi(2,3) = px / (rpred^2);
        Hi_cell{i} = Hi;

        Si(:,:,i) = Hi * Ppred * Hi' + R_polar;

        detSi = det(Si(:,:,i));
        if detSi <= 0
            Si(:,:,i) = Si(:,:,i) + 1e-6 * eye(2);
            detSi = det(Si(:,:,i));
        end
        invSi = inv(Si(:,:,i));
        exponent = -0.5 * (nu(:,i)' * invSi * nu(:,i));
        coef = 1 / (2*pi*sqrt(detSi));
        L(i) = coef * exp(exponent);
    end

    % Validation gating
    validIdx = [];
    for i=1:Nmeas
        md2 = nu(:,i)' * (Si(:,:,i)\nu(:,i));
        if md2 <= gate_thresh
            validIdx(end+1) = i; %#ok<AGROW>
        end
    end

    if isempty(validIdx)
        x = xpred;
        P = Ppred;
        Xest_hist(:,k) = x;
        P_hist(:,:,k) = P;
        assoc_info.beta{k} = [];
        assoc_info.numMeas(k) = 0;
        assoc_info.likelihoods{k} = [];
        continue;
    end

    nu_v = nu(:,validIdx);
    Si_v = Si(:,:,validIdx);
    Hi_v = Hi_cell(validIdx);
    L_v = L(validIdx);
    m = numel(validIdx);

    % Compute PDA association probabilities beta
    if use_lambda_FA
        denom = PD * sum(L_v) + lambda_FA * V_gate;
    elseif ~isempty(lambda_c)
        S_avg = mean(reshape(Si_v,2,2,m),3);
        V_est = pi * gate_thresh * sqrt(det(S_avg));
        denom = PD * sum(L_v) + lambda_c * V_est;
    elseif isfield(params,'lambda_FA') && ~isempty(params.lambda_FA)
        S_avg = mean(reshape(Si_v,2,2,m),3);
        V_est = pi * gate_thresh * sqrt(det(S_avg));
        denom = PD * sum(L_v) + params.lambda_FA * V_est;
    else
        denom = PD * sum(L_v) + 1e-6;
    end

    beta = zeros(m,1);
    for i=1:m
        beta(i) = PD * L_v(i) / denom;
    end
    beta0 = 1 - sum(beta);

    assoc_info.beta{k} = beta;
    assoc_info.numMeas(k) = m;
    assoc_info.likelihoods{k} = L_v;

    % Compute Kalman gains Ki for each validated measurement
    Ki_cell = cell(1,m);
    for i=1:m
        Hi = Hi_v{i};
        Si_i = Si_v(:,:,i);
        Ki_cell{i} = (Ppred * Hi') / Si_i; % 4x2
    end

    % State update: x = xpred + sum_i beta_i * K_i * nu_i
    dx = zeros(4,1);
    for i=1:m
        dx = dx + beta(i) * (Ki_cell{i} * nu_v(:,i));
    end
    x_upd = xpred + dx;

    % Covariance update (practical PDA approximation)
    P_upd = beta0 * Ppred;
    P_spread = zeros(4,4);
    for i=1:m
        Ki = Ki_cell{i};
        Hi = Hi_v{i};
        Ri = R_polar;
        Pi = (I4 - Ki*Hi) * Ppred * (I4 - Ki*Hi)' + Ki * Ri * Ki';
        P_upd = P_upd + beta(i) * Pi;

        x_i_hat = xpred + Ki * nu_v(:,i);
        diff = x_i_hat - x_upd;
        P_spread = P_spread + beta(i) * (diff * diff');
    end
    P_upd = P_upd + P_spread;

    P_upd = (P_upd + P_upd')/2;
    [Vd, Dd] = eig(P_upd);
    Dd = diag(max(diag(Dd), 1e-8));
    P_upd = Vd * diag(Dd) * Vd';

    x = x_upd;
    P = P_upd;

    Xest_hist(:,k) = x;
    P_hist(:,:,k) = P;
end

end
