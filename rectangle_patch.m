function [handle_extent]= rectangle_patch(x, y, length,width, theta,color,linestyle,size,...
    Marker_temp,Markersize_FontSize)
% Function: rectangle_patch
% This function creates a rotated rectangular patch in a MATLAB figure.

rot = [cos(theta),-sin(theta);
    sin(theta),cos(theta)];

xs = [-length/2, -length/2, +length/2, +length/2];
ys = [+width/2, -width/2, -width/2, +width/2];
xy_rot = rot*[xs;ys];
xy_rot = xy_rot + [x;y];
handle_extent = patch('XData',xy_rot(1,:) ,'YData', xy_rot(2,:),'FaceColor','none','EdgeColor',...
    color,'LineWidth',size,'LineStyle',linestyle,'Marker' ...
        ,Marker_temp,'Markersize',Markersize_FontSize);
end

