function [distance] = func_d_gaussian_wasserstein(EO1, EO2)
% Function: func_d_gaussian_wasserstein
% This function calculates the Gaussian Wasserstein distance between two extended objects.

x1 = EO1(1:2);  
theta1 = EO1(3);
eigen_val1 = [EO1(4), EO1(5)];
eigen_vec1 = [cos(theta1), -sin(theta1); sin(theta1), cos(theta1)];
sigma1 = eigen_vec1*diag(eigen_val1.^2)*eigen_vec1';
sigma1 = (sigma1 + sigma1')/2; % make covariance symmetric


x2 = EO2(1:2);  
theta2 = EO2(3);
eigen_val2 = [EO2(4),EO2(5)];
eigen_vec2 = [cos(theta2), -sin(theta2); sin(theta2), cos(theta2)];
sigma2 = eigen_vec2*diag(eigen_val2.^2)*eigen_vec2';
sigma2 = (sigma2 + sigma2')/2; % make covariance symmetric

error_sq = norm(x1-x2)^2 + trace(sigma1 + sigma2 -2*sqrtm((sqrtm(sigma1)*sigma2*sqrtm(sigma1))));
distance = sqrt(error_sq);

end