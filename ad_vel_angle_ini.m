function [xk_1,Pk_1,alphak_1,betak_1,phik_1,phi_sigmak_1,dirkL_1,Q_x,Q_omega,ffactor_IG,ffactor_Dir,...
    C_alphak_1,C_betak_1,C_ffactor_IG] = ad_vel_angle_ini(Para)
% Function: ad_vel_angle_ini
% This function is responsible for initializing the states, parameters, and noise covariances. 
% It extracts the necessary values from 
% the input parameter structure 'Para' and assigns them to the corresponding output variables.

%% Initialization
% xk
xk_1 = Para.ad_vel_angle.xk_1; 
Pk_1 = Para.ad_vel_angle.Pk_1;

% Xk
alphak_1 = Para.ad_vel_angle.alphak_1;
betak_1 = Para.ad_vel_angle.betak_1 ;

% SSF
C_alphak_1 = Para.ad_vel_angle.C_alphak_1 ;
C_betak_1 = Para.ad_vel_angle.C_betak_1 ;


% phik
phik_1 = Para.ad_vel_angle.phik_1 ;
phi_sigmak_1 = Para.ad_vel_angle.phi_sigmak_1;

% pik
dirkL_1 = Para.ad_vel_angle.dirkL_1;

% Q_x
Q_x = Para.ad_vel_angle.Q_x;

% Q_omega
Q_omega = Para.ad_vel_angle.Q_omega;

% forgetting factors
ffactor_IG = Para.ad_vel_angle.ffactor_IG;
ffactor_Dir = Para.ad_vel_angle.ffactor_Dir;
C_ffactor_IG  = Para.ad_vel_angle.C_ffactor_IG;

end