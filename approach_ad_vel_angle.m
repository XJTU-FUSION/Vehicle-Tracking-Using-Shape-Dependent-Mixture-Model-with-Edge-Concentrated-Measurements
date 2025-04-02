function [xk_1,Pk_1,phik_1,phi_sigmak_1,betak_1,alphak_1,...
    C_alphak_1,C_betak_1,dirkL_1,evaluate_index,VB_par] = approach_ad_vel_angle(xk_1,Pk_1,phik_1,phi_sigmak_1,betak_1,alphak_1,...
    C_alphak_1,C_betak_1,dirkL_1,Q_x,Q_omega,ffactor_IG,ffactor_Dir,scale,...
    T,C_ffactor_IG,have_Ns,Ns_2D_all,gt,m,t,y,y_velocity,evaluate_index)
% Function: approach_ad_vel_angle  
% This function executes one cycle of the proposed EOT-SD-V algorithm, integrating prediction, 
% update, and evaluation metric computations.  

H = [1,0,0,0;
    0,1,0,0];
[L_method,h1t,h2t,Ch1t,Ch2t] = set_ad_angle(scale);
if t == 1
    dirkL_1 = 1*ones(1,L_method);
end


%% Prediction
if t == 1
    xkk_1 = xk_1;
    Pkk_1 = Pk_1;
    phikk_1 = phik_1;
    phi_sigmakk_1 = phi_sigmak_1;
    betakk_1 = betak_1;
    alphakk_1 = alphak_1;
    C_betakk_1 = C_betak_1;
    C_alphakk_1 = C_alphak_1;
    dirkkL_1 = dirkL_1;
else
    [xkk_1,Pkk_1,betakk_1,alphakk_1,dirkkL_1] = ...
        state_trans_CT(xk_1,Pk_1,betak_1,alphak_1,dirkL_1,ffactor_IG,ffactor_Dir,...
        Q_x,L_method,phik_1,T);

    F_phi = [1,T;0,1];
    G_phi = [T;1];
    phikk_1 = F_phi*phik_1;
    phi_sigmakk_1 = F_phi*phi_sigmak_1*F_phi'+G_phi*Q_omega*G_phi';

    C_betakk_1 = C_ffactor_IG*C_betak_1;
    C_alphakk_1 = C_ffactor_IG*C_alphak_1;
end

R_a = (0.15*pi/180)^2;
R_r = (0.02)^2;
R_part1 = diag([R_r,R_a]);
position_k_global = H*xkk_1;
if have_Ns
    Ns_temp = Ns_2D_all(:,t);
    position_k = position_k_global-Ns_temp;
else
    position_k = position_k_global;
end
rk_pos = (position_k(1)^2+position_k(2)^2)^(0.5);
cos_pos = position_k(1)/rk_pos;
sin_pos = position_k(2)/rk_pos;
R_part2 = [cos_pos,-position_k(2);
    sin_pos,position_k(1)];
R = R_part2*R_part1*R_part2';

%% Update
var_velocity = 4;
[xk_o,Pk_o,betak_o,alphak_o,phik_o,phi_sigmak_o,dirkL_o,C_betak_o,C_alphak_o] = ...
    Update_ad_vel_Angle(xkk_1,Pkk_1,betakk_1,alphakk_1,C_betakk_1,C_alphakk_1,phikk_1,...
    phi_sigmakk_1,dirkkL_1,y,H,R,L_method,h1t,h2t,...
    y_velocity,Ns_temp,var_velocity,Ch1t,Ch2t);
% xk
xk_1 = xk_o;    Pk_1 = Pk_o;

% Xk
betak_1 = betak_o; alphak_1 = alphak_o;

% pik
dirkL_1 = dirkL_o;

% SSF
C_betak_1 = C_betak_o;
C_alphak_1 = C_alphak_o;

% len_VB
len_VB = (betak_1./(alphak_1-1)).^0.5;
evaluate_index.ad_vel_angle.len_est_all{m,t} = len_VB;


% thetak
phik_1 = phik_o;
phi_sigmak_1 = phi_sigmak_o;
ang_VB = phik_1(1,1);
evaluate_index.ad_vel_angle.ang_all(m,t) =ang_VB;
evaluate_index.ad_vel_angle.sigma_ang_all(m,t) = phi_sigmak_1(1,1);
VB_par = [H*xk_1; len_VB; ang_VB];
evaluate_index.ad_vel_angle.VB_estimate_algorithm(:,t,m) = VB_par;

%% Evaluation metric computations
gt_frame = gt(:, t);

% GWD
ellipse_estimate = [xk_o(1),xk_o(2),ang_VB,len_VB(1),len_VB(2)];
ellipse_truth = [gt_frame(1),gt_frame(2),gt_frame(3),gt_frame(4),gt_frame(5)];
[gw_distance_frame] = func_d_gaussian_wasserstein(ellipse_truth, ellipse_estimate);
evaluate_index.ad_vel_angle.gw_distance_all(m,t) = gw_distance_frame;

% position SE
Position_SE = (xk_o(1)-gt_frame(1))^2+(xk_o(2)-gt_frame(2))^2;
Position_RSE = Position_SE.^(0.5);
evaluate_index.ad_vel_angle.Position_SE_all(m,t) = Position_SE;
evaluate_index.ad_vel_angle.Position_RSE_all(m,t) = Position_RSE;

% orientation SE
abs_error = abs(ang_VB-gt_frame(3));
angle_SE = (abs_error)^2;
angle_RSE = abs_error;
evaluate_index.ad_vel_angle.angle_SE_all(m,t) =  angle_SE;
evaluate_index.ad_vel_angle.angle_RSE_all(m,t) = angle_RSE;

end