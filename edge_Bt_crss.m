function [Bt_edge, edgecrss_output] = edge_Bt_crss(T,a_alloy,flag_plot)

edgeburger = a_alloy/2*sqrt(3);% A


%% edge friction_coefficient
edgestress_300k = [0.5036,0.6055,0.7036,0.8073,0.9091,1.0109];
edgevelocity_300k = [4.4444,5.3968,6.6984,8.2857,9.4601,10.7619];

edgestress_1500k = [0.5036,0.6055,0.7055,0.8054,0.9091,1.0091];
edgevelocity_1500k = [3.873,5.3651,6.6032,8.3492,9.8094,10.7936];

edgestress_2000k = [0.5036,0.6036,0.7055,0.8073,0.9073,1.0091];
edgevelocity_2000k = [3.9048,4.6032,6.4762,7.8094,8.9841,10.0952];


stress_tensor = [edgestress_300k;edgestress_1500k;edgestress_2000k];
velocity_tensor = [edgevelocity_300k;edgevelocity_1500k;edgevelocity_2000k];



%% calculate friction coefficient
edgestress_avg = (sum(stress_tensor)/size(stress_tensor,1));  % GPa
edgevelocity_avg = (sum(velocity_tensor)/size(stress_tensor,1));  %A/ps

coefficients = polyfit(edgestress_avg,edgevelocity_avg, 1);
edge_Bt_fitted = 1/coefficients(1);
edge_displacement = coefficients(2);
edge_velocity_fitted = edgestress_avg/edge_Bt_fitted+edge_displacement;

Bt_edge_GPa = edge_Bt_fitted*edgeburger*1e-12;  %GPa*s
Bt_edge = Bt_edge_GPa*1e-11; %N/A^2*s

%% CRSS fitting
edgethreshold_T = [300,600,800,1500,1800,2000];
edgecrss = [0.3509,0.3509,0.3509,0.2509,0.2509,0.1];
edgecrss_coefficient = polyfit(edgethreshold_T,edgecrss, 1);
edgecrss_fitted = edgecrss_coefficient(1)*edgethreshold_T+edgecrss_coefficient(2);
    
edgecrss_GPa = edgecrss_coefficient(1)*T+edgecrss_coefficient(2);    % in GPa

edgecrss_output = edgecrss_GPa*1e-11; % N/A^2

end
