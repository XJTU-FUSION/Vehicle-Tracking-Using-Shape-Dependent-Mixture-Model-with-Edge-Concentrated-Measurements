function [xk_vel_cov,xk_vel_mean] =...
    VB_vel_xk(EZZ,G1,G2,GZ2,GZ4,zk,y_velocity,var_vel,Ms,Nmeas,omegak,omega_sigmak)
% Function: VB_vel_xk  
% This function can obtain the contribution of Doppler velocity
% measurements when calculating q_x

% covariance
xk_vel_cov = 0;
for jj = 1:Nmeas
    EZZ_temp = EZZ{1,jj};
    temp1 = G1*EZZ_temp*(G1)';
    temp2 = omegak*G1*EZZ_temp*G2;
    temp3 = omegak*G2'*EZZ_temp*G1';
    temp4 = ((omegak)^2+omega_sigmak)*G2'*EZZ_temp*G2;
    temp_cov = (1/var_vel)*(temp1+temp2+temp3+temp4); 
    xk_vel_cov = xk_vel_cov + temp_cov;
end

% mean
xk_vel_mean = 0;
for jj = 1:Nmeas
    EZZ_temp = EZZ{1,jj};
    GZ2_temp = GZ2{1,jj};
    GZ4_temp = GZ4{1,jj};
    part1 = omegak*Ms'*EZZ_temp*G1'; 
    part2 = ((omegak)^2+omega_sigmak)*Ms'*EZZ_temp*G2; 
    mix1 = (-1/var_vel)*((part1+part2)'); 


    part3 = G1*GZ4_temp*zk(:,jj);
    part4 = omegak*G2'*GZ4_temp*zk(:,jj);
    part5 = G1*GZ2_temp'; 
    part6 = omegak*G2'*GZ2_temp';
    mix2= (y_velocity(jj)/(var_vel))*(part3+part4+part5+part6);

    temp_mean =mix1+mix2;
    xk_vel_mean = xk_vel_mean + temp_mean;
end
end