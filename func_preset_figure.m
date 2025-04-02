function [Marker_all,Marker_temp_gt,LineStyle_gt,Linewidth_FontSize,Linewidth_gt_FontSize,Markersize_FontSize,xylabel_Fontsize,legend_FontSize,...
    estimate_ellip_size,mea_size,color_mea,color_gt,color_all,legend_cell,legend_cell_RMSE,linestyle_all,algorithm_mode_all] = func_preset_figure
% Function: func_preset_figure
% This function is used to preset various parameters for plotting figures, such as marker styles, line styles, 
% colors, sizes, and legend labels. These parameters can be used to customize the appearance of plots, 
% including squares, measurements, and ground truth data.

Marker_all = {};
Marker_temp_gt = '+';
LineStyle_gt = '-';
%% 
Linewidth_FontSize = 1.5;
Linewidth_gt_FontSize = 2;
Markersize_FontSize = 6;
xylabel_Fontsize = 16;
legend_FontSize = 14;
%% 
estimate_ellip_size = 2;
mea_size = 4;
%%
color_mea = 'b*';
color_gt = [99/255,227/255,152/255]; 
color_all = {};
legend_cell = {'Ground truth','Measurement'};
legend_cell_RMSE = cell(1,1);
linestyle_all = cell(1,1);
algorithm_mode_all = cell(1,1);

end