function [xkk_1,Pkk_1,betakk_1,alphakk_1,dirkkL_1] = ...
    state_trans_CT(xk_1,Pk_1,betak_1,alphak_1,dirkL_1,ffactor_IG,ffactor_Dir,...
    Q_x,L_method,phik_1,T)
% Function: state_trans_CT  
% This function performs the prediction step of the EOT-SD-V algorithm.


omega_e = phik_1(2,1);
F_x_CT = zeros(4,4);
sin_temp = sin(omega_e*T);
cos_temp = cos(omega_e*T);
sin_temp2 = sin_temp/omega_e;
cos_temp2 = (1-cos_temp)/omega_e;
F_x_CT(1,:) = [1,0,sin_temp2,-cos_temp2];
F_x_CT(2,:) = [0,1,cos_temp2,sin_temp2];
F_x_CT(3,:) = [0,0,cos_temp,-sin_temp];
F_x_CT(4,:) = [0,0,sin_temp,cos_temp];

% xk
xkk_1 = F_x_CT*xk_1;
Pkk_1 = F_x_CT*Pk_1*F_x_CT' + Q_x;

% Xk
betakk_1 = ffactor_IG*betak_1;
alphakk_1 = ffactor_IG*alphak_1;

% pik
dirkkL_1 = zeros(1,L_method);
for l = 1:L_method
    dirkkL_1(1,l) = ffactor_Dir*dirkL_1(1,l);
end
end