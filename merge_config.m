clear;clc;close all;

% load 1jog.mat
% subplot(2,2,1)
% plotnodes_alan_normalizedscale(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,typelist,viewangle,plim_x,plim_y,plim_z,includeSuperjogs); 
% 
% 
% load 2jog.mat
% hold on 
% subplot(2,2,2)
% 
% plotnodes_alan_normalizedscale(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,typelist,viewangle,plim_x,plim_y,plim_z,includeSuperjogs); 
% hold off
% 
% load 3jog.mat
% subplot(2,2,3)
% 
% hold on 
% plotnodes_alan_normalizedscale(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,typelist,viewangle,plim_x,plim_y,plim_z,includeSuperjogs); 
% hold off

load output_1000stress_2000T_1636step.mat
% subplot(2,2,4)
plim_y=100
hold on 

plotnodes_alan_normalizedscale(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,typelist,viewangle,plim_x,plim_y,plim_z,includeSuperjogs); 
hold off

% % 
% load 0jog_1700.mat
% plim_y=100
% 
% hold on 
% plotnodes_alan_normalizedscale(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,typelist,viewangle,plim_x,plim_y,plim_z,includeSuperjogs); 
% hold off



