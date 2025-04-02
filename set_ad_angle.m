function [L_method,h1t_all_out,h2t_all_out,Ch1t_all_out,Ch2t_all_out] = set_ad_angle(scale)
% Function: set_ad_angle
% This function is used to calculate h_t in (19) and Cht in (20).
h1t_1 = [1,0]';
h1t_2 = [0,1]';
h1t_3 = [-1,0]';
h1t_4 = [0,-1]';

h1t_sample = [h1t_1,h1t_2,h1t_3,h1t_4];
h1t_all = h1t_sample;

h2t_temp = -h1t_sample;
h2t_all = h2t_temp;


%% Cht
Ch1t_1 = scale*diag([1,0]);
Ch1t_3 = Ch1t_1;
Ch1t_2 = scale*diag([0,1]);
Ch1t_4 = Ch1t_2;

Ch2t_1 = scale*diag([0,1]);
Ch2t_3 = Ch2t_1;
Ch2t_2 = scale*diag([1,0]);
Ch2t_4 = Ch2t_2;

L_method = 4;
h1t_all_out = h1t_all;
h2t_all_out = h2t_all;
Ch2t_all_out = {Ch2t_1,Ch2t_2,Ch2t_3,Ch2t_4};
Ch1t_all_out ={Ch1t_1,Ch1t_2,Ch1t_3,Ch1t_4};


end


