function X = generateTarget_CV(x0, K, dt, Qproc_xy)
% 生成真实目标轨迹（CV 模型）
A = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1];
X = zeros(4,K);
X(:,1) = x0;
for k=2:K
    eps = mvnrnd([0 0], Qproc_xy)'; % [epsx; epsy]
    proc = [0.5*eps(1)*dt^2; eps(1)*dt; 0.5*eps(2)*dt^2; eps(2)*dt];
    X(:,k) = A*X(:,k-1) + proc;
end
end
