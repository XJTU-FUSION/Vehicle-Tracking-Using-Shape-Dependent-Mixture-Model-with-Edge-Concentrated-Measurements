%% Vehicle Tracking Using Shape-Dependent Mixture Model with Edge-Concentrated Measurements
% reference:
% "Vehicle Tracking Using Shape-Dependent Mixture Model with Edge-Concentrated Measurements"
% (IEEE Transactions on Intelligent Transportation Systems) by Zheng Wen, Jian Lan*, Le Zheng, Tao Zeng, 2025.
% Copyright by Zheng Wen, Beijing Institute of Technology; Jian Lan, Xi'an Jiaotong University
% Email:  1459223362@qq.com; lanjian@mail.xjtu.edu.cn

close all
clc
clear
dbstop warning
% preset parameters for plotting figures
[Marker_all,Marker_temp_gt,LineStyle_gt,Linewidth_FontSize,Linewidth_gt_FontSize,Markersize_FontSize,xylabel_Fontsize,legend_FontSize,...
    estimate_ellip_size,mea_size,color_mea,color_gt,color_all,legend_cell,legend_cell_RMSE,linestyle_all,algorithm_mode_all] = func_preset_figure;
set(gca,'looseInset',[0.01,0.01,0.01,0.01])

Mc_num = 1;
visualize = 1;
frame_view = 1;
evaluate_index = [];
exist_gt = 0;
algorithm_num = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Started
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for algorithm = 1:algorithm_num
    %% Algorithm information
    mode= 'ad_vel_angle';
    legend_name ='EOT-SD-V';
    color_now = [255/255,0/255,0/255];
    Marker_now = '*';
    linestyle_temp = '-';

    scale = 1/3;
    T = 0.1; % time interval
    frame_num = 42; % frame number
    have_Ns = 1;
    Ns_2D_all = load("Ns_demo.mat").Ns_2D_all; % N_s
    gt = load('gt_demo').gt; % ground truth
    Para_set = load('Para_demo.mat').Para_set; % parameter setting for the algorithm 

    for m = 1:Mc_num
        %% Initialization
        [xk_1,Pk_1,alphak_1,betak_1,phik_1,phi_sigmak_1,...
            dirkL_1,Q_x,Q_omega,ffactor_IG,ffactor_Dir,C_alphak_1,C_betak_1,C_ffactor_IG] = ad_vel_angle_ini(Para_set);
        Measurement = load('mea_data_demo.mat').mea_data;
        for t = 1:frame_num
            %% Load measurement for each frame
            data_frame = Measurement{1,t};
            data_x = data_frame(1,:);
            data_y = data_frame(2,:);
            y = [data_x;data_y];
            y_velocity = data_frame(3,:);

            %% Filtering
            [xk_1,Pk_1,phik_1,phi_sigmak_1,betak_1,alphak_1,...
                C_alphak_1,C_betak_1,dirkL_1,evaluate_index,VB_par] = approach_ad_vel_angle(xk_1,Pk_1,phik_1,phi_sigmak_1,betak_1,alphak_1,...
                C_alphak_1,C_betak_1,dirkL_1,Q_x,Q_omega,ffactor_IG,ffactor_Dir,scale,...
                T,C_ffactor_IG,have_Ns,Ns_2D_all,gt,m,t,y,y_velocity,evaluate_index);
            %% Plot estimate
            color_all{1,algorithm} = color_now;
            Marker_all{1,algorithm} = Marker_now;
            Marker_temp = Marker_now;
            color_temp = color_all{1,algorithm};
            if visualize
                if mod(t,frame_view)==0
                    t
                    hold on
                    axis equal
                    meas_points = plot(y(1,:), y(2,:), color_mea,'LineWidth',mea_size,'MarkerSize',1);
                    gt_size = 2;
                    if ~exist_gt
                        gt_plot = rectangle_patch(gt(1, t), gt(2, t), 2*gt(4, t),2*gt(5, t), gt(3, t),color_gt,'-',gt_size,...
                            Marker_temp_gt,Markersize_FontSize);
                    end
                    method1 = rectangle_patch(VB_par(1), VB_par(2), 2*VB_par(3),2*VB_par(4), VB_par(5),color_temp,linestyle_temp,gt_size,...
                        Marker_temp,Markersize_FontSize);
                    pause(0.1);
                    xlabel('x (m)','fontsize',14);
                    ylabel('y (m)','fontsize',14);
                end
            end
        end
    end
    exist_gt = 1;
    legend_cell{1,2+algorithm} = legend_name;
    legend_cell_RMSE{1,algorithm} = legend_name;
    linestyle_all{1,algorithm} = linestyle_temp;
    algorithm_mode_all{1,algorithm} = mode;
end
if visualize
    box on
    box(gca,'on');
    legend([gt_plot,meas_points,method1], legend_cell,'Location','northwest','fontsize',12,'NumColumns',1);
end
hold off
set(gca,'Xcolor',[0 0 0]);
set(gca,'Ycolor',[0 0 0]);
marker_idx = 1:5:frame_num;

%% Showing results
% GWD
figure()
set(gca,'looseInset',[0.01,0.01,0.01,0.01])
for algorithm =1:algorithm_num
    mode_name = algorithm_mode_all{1,algorithm};
    gw_distance_all = evaluate_index.(mode_name).gw_distance_all;
    [gw_distance_RMSE_all] = func_cal_RMSE(gw_distance_all);
    color_temp = color_all{1,algorithm};
    linestyle_temp = linestyle_all{1,algorithm};
    Marker_temp = Marker_all{1,algorithm};
    plt_pos = plot(1:frame_num,gw_distance_RMSE_all(1:frame_num),'LineWidth',Linewidth_FontSize,...
        'Color',color_temp,'LineStyle',linestyle_temp,'Marker' ...
        ,Marker_temp,'Markersize',Markersize_FontSize,'MarkerIndices',marker_idx);

    hold on
    xlabel('Time Step','fontsize',xylabel_Fontsize)
    ylabel('Average GWD (m)','fontsize',xylabel_Fontsize)
end
xlim([1,frame_num]);
ylim([0,1.8]);
set(legend,'Location','NorthWest')
legend(legend_cell_RMSE,'fontsize',legend_FontSize,'NumColumns',2);

% Position RMSE
figure()
set(gca,'looseInset',[0.01,0.01,0.01,0.01])
for algorithm =1:algorithm_num
    mode_name = algorithm_mode_all{1,algorithm};
    Position_RSE_all = evaluate_index.(mode_name).Position_RSE_all;
    [Position_RMSE_all] = func_cal_RMSE(Position_RSE_all);
    color_temp = color_all{1,algorithm};
    linestyle_temp = linestyle_all{1,algorithm};
    Marker_temp = Marker_all{1,algorithm};
    plt_pos = plot(1:frame_num,Position_RMSE_all(1:frame_num),'LineWidth',Linewidth_FontSize,'Color',color_temp,'LineStyle',linestyle_temp,'Marker' ...
        ,Marker_temp,'Markersize',Markersize_FontSize,'MarkerIndices',marker_idx);

    hold on
    xlabel('Time Step','fontsize',xylabel_Fontsize)
    ylabel('Centroid Position RMSE (m)','fontsize',xylabel_Fontsize)
end
xlim([1,frame_num]);
ylim([0,1.7]);
set(legend,'Location','NorthWest')
legend(legend_cell_RMSE,'fontsize',legend_FontSize,'NumColumns',2);

% Orientation RMSE
figure()
set(gca,'looseInset',[0.01,0.01,0.01,0.01])
for algorithm =1:algorithm_num
    mode_name = algorithm_mode_all{1,algorithm};
    angle_RSE_all = evaluate_index.(mode_name).angle_RSE_all;
    [angle_RMSE_all] = func_cal_RMSE(angle_RSE_all);
    color_temp = color_all{1,algorithm};
    linestyle_temp = linestyle_all{1,algorithm};
    Marker_temp = Marker_all{1,algorithm};
    plt_angle = plot(1:frame_num,angle_RMSE_all(1:frame_num),'LineWidth',Linewidth_FontSize,'Color',color_temp,'LineStyle',linestyle_temp,'Marker' ...
        ,Marker_temp,'Markersize',Markersize_FontSize,'MarkerIndices',marker_idx);

    hold on
    xlabel('Time Step','fontsize',xylabel_Fontsize)
    ylabel('Orientation RMSE (rad)','fontsize',xylabel_Fontsize)
end
xlim([1,frame_num]);
ylim([0,0.65]);
set(legend,'Location','NorthWest')
legend(legend_cell_RMSE,'fontsize',legend_FontSize,'NumColumns',2);





