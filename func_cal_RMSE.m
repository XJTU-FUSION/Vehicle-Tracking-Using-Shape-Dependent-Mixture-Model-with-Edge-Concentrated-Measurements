function [Position_RMSE_all] = func_cal_RMSE(Position_RSE_all)
% Function: func_cal_RMSE
% This function calculates the root mean squared error (RMSE) of the centroid position.

Position_SE_all = Position_RSE_all.^2;
Position_MSE_all = mean(Position_SE_all,1);
Position_RMSE_all = Position_MSE_all.^0.5;

end