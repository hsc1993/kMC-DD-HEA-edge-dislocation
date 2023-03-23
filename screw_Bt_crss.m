function [Bt_screw, screwcrss_output] = screw_Bt_crss(T,a_alloy,flag_plot)

screwburger = a_alloy/2*sqrt(3);% A

%% screw friction_coefficient
screw_temperature = [300;600;800;1000;1200;1800;2000];
screwstress_300k = [1.59,1.795,2.19,2.39,2.59,2.99];
screwvelocity_300k = [0.089,0.615,1.96,2.476,2.91,3.65];

screwstress_600k = [1.195,1.395,1.59,2.19,2.586,2.99];
screwvelocity_600k = [0.124,0.476,0.741,2.037,3.024,3.554];

screwstress_800k = [1.395,1.79,1.995,2.19,2.59,2.99];
screwvelocity_800k = [0.56,1.65,2.27,2.66,2.63,3.75];


screwstress_1000k = [0.995,1.60,1.79,2.19,2.59,2.99];
screwvelocity_1000k = [0.391,1.513,1.66,2.31,3.187,3.87];

screwstress_1200k = [0.995,1.195,1.595,2,2.59,2.79];
screwvelocity_1200k = [0.2825,0.776,1.34,2.17,3.31,3.724];

screwstress_1800k = [0.6,0.995,1.2,1.395,1.60,2.39];
screwvelocity_1800k = [0.2522,0.5,0.84,0.846,1.44,2.98];

screwstress_2000k = [1.795,1.995,2.19,2.39,2.59,2.795];
screwvelocity_2000k = [2.3,2.533,2.809,2.99,3.437,3.985];


stress_tensor = [screwstress_300k;screwstress_600k;screwstress_800k;screwstress_1000k;screwstress_1200k;screwstress_1800k;screwstress_2000k];
velocity_tensor = [screwvelocity_300k;screwvelocity_600k;screwvelocity_800k;screwvelocity_1000k;screwvelocity_1200k;screwvelocity_1800k;screwvelocity_2000k];


%% calculate friction coefficient



screw_Bt_fitted = zeros(size(stress_tensor,1),1);
screw_displacement = zeros(size(stress_tensor,1),1);
screw_Bt_gpa = zeros(size(stress_tensor,1),1);

for ii = 1:size(stress_tensor,1)
    coefficients = polyfit(stress_tensor(ii,:),velocity_tensor(ii,:), 1);
    screw_Bt_fitted(ii) = 1/coefficients(1);  %GPa*ps/A
    screw_displacement(ii) = coefficients(2);
    screw_Bt_gpa(ii) = screw_Bt_fitted(ii)*screwburger*1e-12;  %GPa*s
end


fun = @(x,xdata) x(1)*log(x(2)*xdata)+x(3);
x0 = [1,1,0];
screw_Bt_log_fit_params = lsqcurvefit(fun,x0,screw_temperature,screw_Bt_gpa);
screw_Bt_log_fitted = fun(screw_Bt_log_fit_params,screw_temperature);
Bt_screw_GPa = fun(screw_Bt_log_fit_params,T);
Bt_screw = Bt_screw_GPa*1e-11; %N/A^2*s


%% CRSS fitting
screwthreshold_T = [300,600,800,1000,1200,1800,2000];
screwcrss = [1.395,0.995,0.795,0.405,0.395,0.195,0];
screwcrss_coefficient = polyfit(screwthreshold_T,screwcrss, 1);
screwcrss_fitted = screwcrss_coefficient(1)*screwthreshold_T+screwcrss_coefficient(2);
screwcrss_GPa = screwcrss_coefficient(1)*T+screwcrss_coefficient(2);    % in GPa

screwcrss_output = screwcrss_GPa*1e-11; % N/A^2

end
