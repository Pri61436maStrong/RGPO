function [xi, MD] = computeScore_strict(Xtrue, Xest, assoc_info, gate_thresh)
% 计算 xi (轨迹中断率) 与 MD (平均误差距离)
K = size(Xtrue,2);
posTrue = Xtrue([1,3],:);
posEst = Xest([1,3],:);
dists = sqrt(sum((posTrue - posEst).^2,1)); % per-step position errors
MD = mean(dists); % 平均误差距离

% 航迹中断率 xi:
lost_count = 0;
for k=1:K
    betas = assoc_info.beta{k};
    if isempty(betas)
        lost_count = lost_count + 1;
    else
        if sum(betas) < 0.5
            lost_count = lost_count + 1;
        end
    end
end
xi = lost_count / K; % 轨迹中断率
end
