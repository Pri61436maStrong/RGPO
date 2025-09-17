function Xf = generateFakeFromStrategy(Xr, dtaus, dt)
% Xr: 4xK true states
% dtaus: 1xK array of Δτ_k in seconds
K = size(Xr,2);
tau = cumsum(dtaus); % cumulative delay at each step
c = 3e8;
dk = c * tau / 2;    % radial offset (meters)
Xf = Xr;
for k=1:K
    xr = Xr(1,k); yr = Xr(3,k);
    r = sqrt(xr^2 + yr^2);
    if r < 1e-6
        theta = 0;
    else
        theta = atan2(yr, xr);
    end
    dx = dk(k) * cos(theta);
    dy = dk(k) * sin(theta);
    Xf(1,k) = Xr(1,k) + dx;
    Xf(3,k) = Xr(3,k) + dy;
end
end
