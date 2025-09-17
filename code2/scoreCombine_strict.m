function r = scoreCombine_strict(xi, MD, w1, w2)
% 按论文式(2-29) 线性加权组合 r = w1*xi + w2*MD
r = w1 * xi + w2 * MD;
end
