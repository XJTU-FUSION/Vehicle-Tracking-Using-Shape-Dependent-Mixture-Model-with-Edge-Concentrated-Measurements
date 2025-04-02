function [ln_Cht] = cal_ln_Cht(a,b,scale)
% Function: cal_ln_Cht  
% This function calculates the expectation of ln_Cht.

t =1;
p1 = gammainc(b/t,a,'upper'); 
fun = @(x) (gammainc(x,a,'upper'))./p1./x;
ln_Cht =log((scale^2)*b)-log(b/t) -integral(fun,b/t,Inf);


end