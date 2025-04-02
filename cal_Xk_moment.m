function [s_estimate_out] = cal_Xk_moment(a_all,b_all,k)
% Function: cal_Xk_moment
% This function calculates the k-th moment of Xk.

s_estimate_out = zeros(2,2);
for n =1:2
    a = a_all(n,1);
    b = b_all(n,1);
    log_1 = gammaln(a-k);
    log_2 = gammaln(a);

    log_mix = log_1-log_2;
    part1 = exp(log_mix);
    part2 = b^k;
    s_estimate_out(n,n) = part1*part2;
end