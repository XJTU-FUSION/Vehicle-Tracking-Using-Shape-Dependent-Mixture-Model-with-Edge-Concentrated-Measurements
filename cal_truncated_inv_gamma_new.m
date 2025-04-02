function [s_estimate_out] = cal_truncated_inv_gamma_new(a,b,k)
% Function: cal_truncated_inv_gamma_new
% This function calculates the k-th moment of Ck.

t=1;
a1 = a(1);a2 = a(2);
b1 = b(1);b2 = b(2);

if k>0
    temp_new1 = beta(a1-k,k)/gamma(k);
    temp_new2 = beta(a2-k,k)/gamma(k);
else
    temp_new1 = gamma(-k)/beta(a1,-k);
    temp_new2 = gamma(-k)/beta(a2,-k);
end

value11= gammainc(b1/t,a1-k,'upper');
value12 = gammainc(b1/t,a1,'upper');

value21= gammainc(b2/t,a2-k,'upper');
value22 = gammainc(b2/t,a2,'upper');

s_estimate1 = temp_new1*(b1^(k))*(value11/value12);
s_estimate2 = temp_new2*(b2^(k))*(value21/value22);

s_estimate_out = [s_estimate1,0;0,s_estimate2];

end