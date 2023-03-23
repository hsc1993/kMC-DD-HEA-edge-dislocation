% function [a_alloy,Bt_edge_GPa,Bt_screw_GPa,debye_frequency,edgecrss,screwcrss] = hea_constants(T,flag_plot)
clc;clear all;close all
flag_plot=1;
T=300;

% physical constants
kb = 8.617333262145e-5;	% eV * K^-1
h_bar = 6.582119569e-16; % eV * s

% lattice constant in A
a_Nb = 3.30;
a_Mo = 3.15;
a_Ta = 3.31;
a_W = 3.16;
a_alloy = 0.25*(a_Nb+a_Mo+a_Ta+a_W);
edgeburger = a_alloy/2*sqrt(3);% A
screwburger = a_alloy/2*sqrt(3);% A

%% debye_frequency_alloy
debye_T_Nb = 275;
debye_T_Mo = 450;
debye_T_Ta = 240;
debye_T_W = 400;

debye_frequency_Nb = kb*debye_T_Nb/h_bar;  % s^-1
debye_frequency_Mo = kb*debye_T_Mo/h_bar;
debye_frequency_Ta = kb*debye_T_Ta/h_bar;
debye_frequency_W = kb*debye_T_W/h_bar;
debye_frequency_alloy = 0.25*(debye_frequency_Nb+debye_frequency_Mo+debye_frequency_Ta+debye_frequency_W);



%% edge friction_coefficient
edgestress_300k = [0.5036,0.6055,0.7036,0.8073,0.9091,1.0109];
edgevelocity_300k = [4.4444,5.3968,6.6984,8.2857,9.4601,10.7619];

edgestress_1500k = [0.5036,0.6055,0.7055,0.8054,0.9091,1.0091];
edgevelocity_1500k = [3.873,5.3651,6.6032,8.3492,9.8094,10.7936];

edgestress_2000k = [0.5036,0.6036,0.7055,0.8073,0.9073,1.0091];
edgevelocity_2000k = [3.9048,4.6032,6.4762,7.8094,8.9841,10.0952];


edgestress_tensor = [edgestress_300k;edgestress_1500k;edgestress_2000k];
edgevelocity_tensor = [edgevelocity_300k;edgevelocity_1500k;edgevelocity_2000k];


%% screw friction_coefficient
screw_sampling_temperature = [300,600,800,1000,1200,1800,2000];
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


screwstress_tensor = [screwstress_300k;screwstress_600k;screwstress_800k;screwstress_1000k;screwstress_1200k;screwstress_1800k;screwstress_2000k];
screwvelocity_tensor = [screwvelocity_300k;screwvelocity_600k;screwvelocity_800k;screwvelocity_1000k;screwvelocity_1200k;screwvelocity_1800k;screwvelocity_2000k];

edgerealdatarow = size(edgestress_tensor,1);
emptyrow = size(screwstress_tensor,1)-size(edgestress_tensor,1);

stress_tensor = [edgestress_tensor;zeros(emptyrow,size(edgestress_tensor,2))];
stress_tensor(:,:,2) = screwstress_tensor;

velocity_tensor = [edgevelocity_tensor;zeros(emptyrow,size(edgestress_tensor,2))];
velocity_tensor(:,:,2) = screwvelocity_tensor;

%% calculate friction coefficient
edgestress_avg = (sum(stress_tensor(1:edgerealdatarow,:,1))/size(stress_tensor(1:edgerealdatarow,:,1),1));  % GPa
edgevelocity_avg = (sum(velocity_tensor(1:edgerealdatarow,:,1))/size(stress_tensor(1:edgerealdatarow,:,1),1));  %A/ps

coefficients = polyfit(edgestress_avg,edgevelocity_avg, 1);

edge_coefficient_plot = 1/coefficients(1);
edge_displacement = coefficients(2);

edgefriction_coefficient_gpa = edge_coefficient_plot*edgeburger*1e-12;  %GPa*s



screw_coefficient_plot = zeros(size(stress_tensor(:,:,2),1),1);
screwfriction_coefficient_gpa = zeros(size(stress_tensor(:,:,2),1),1);
screw_displacement = zeros(size(stress_tensor(:,:,2),1),1);

for ii = 1:size(stress_tensor(:,:,2),1)
    coefficients = polyfit(stress_tensor(ii,:,2),velocity_tensor(ii,:,2), 1);
    screw_coefficient_plot(ii) = 1/coefficients(1);  %GPa*ps/A
    screw_displacement(ii) = coefficients(2);

    screwfriction_coefficient_gpa(ii) = screw_coefficient_plot(ii)*screwburger*1e-12;  %GPa*s
end
mobility_fitting_coefficients = polyfit(screw_sampling_temperature',screwfriction_coefficient_gpa, 2);


% screwfriction_coefficient_plot = screwstress_avg/screwvelocity_avg;  %GPa*ps/A
% screwfriction_coefficient_gpa = screwstress_avg/screwvelocity_avg*screwburger*1e-12;  %GPa*s



edgethreshold_T = [300,600,800,1500,1800,2000];
edgecrss = [0.3509,0.3509,0.3509,0.2509,0.2509,0.1];
edgecrss_coefficient = polyfit(edgethreshold_T,edgecrss, 1);

screwthreshold_T = [300,600,800,1000,1200,1800,2000];
screwcrss = [1.395,0.995,0.795,0.405,0.395,0.195,0];
screwcrss_coefficient = polyfit(screwthreshold_T,screwcrss, 1);


if flag_plot == 1

    f = figure;
    f.Position = [100 100 1000 800];
    
    for ii = 1:2
        subplot(2,3,ii)
        hold on

        if ii == 1
    %         ii = 1, edge dislocation
            plot(sum(stress_tensor(1:edgerealdatarow,:,ii))/edgerealdatarow,sum(stress_tensor(1:edgerealdatarow,:,ii))/edgerealdatarow./edge_coefficient_plot+edge_displacement,'LineWidth',5)
        end
        

        if ii == 2
    %         ii = 2, screw dislocation
            for jj = 1:size(stress_tensor(:,:,ii),1)
                plot(stress_tensor(jj,:,ii),stress_tensor(jj,:,ii)./screw_coefficient_plot(jj)+screw_displacement(jj),'LineWidth',5)
            end
        end
        hold off

        
        xlabel('stress [GPa]', 'FontSize', 20)
        ylabel('velocity [A/Ps]', 'FontSize', 20)

        ax=gca;
        ax.FontSize = 20;
        hold on
        for i = 1:size(stress_tensor(:,:,ii),1)
            if ii == 1 && i<=edgerealdatarow
                scatter(stress_tensor(i,:,ii),velocity_tensor(i,:,ii),'r')
            end
            if ii == 2
                scatter(stress_tensor(i,:,ii),velocity_tensor(i,:,ii),'r')
            end
        end
        hold off
        
%         add legends
        if ii ==1
            LG=legend(strcat('friction\_coefficient=',num2str(edgefriction_coefficient_gpa),'GPa*s'));
            LG.FontSize = 14;
        end 
        
        if ii ==2
            lg1 = strcat('friction\_coefficient=',num2str(screwfriction_coefficient_gpa(1)),'GPa*s');
            lg2 = strcat('friction\_coefficient=',num2str(screwfriction_coefficient_gpa(2)),'GPa*s');
            lg3 = strcat('friction\_coefficient=',num2str(screwfriction_coefficient_gpa(3)),'GPa*s');
            LG=legend(lg1,lg2,lg3);
            LG.FontSize = 14;
        end 
        
    end


    %% threshold rss
    for ii = 3:4
        subplot(2,3,ii)
        hold on
        if ii == 3
            plot(edgethreshold_T,edgecrss,'LineWidth',2)
            scatter(edgethreshold_T,edgecrss)
            plot(edgethreshold_T,edgecrss_coefficient(1)*edgethreshold_T+edgecrss_coefficient(2),'LineWidth',2)
        else
            plot(screwthreshold_T,screwcrss,'LineWidth',2)
            scatter(screwthreshold_T,screwcrss)  
            plot(screwthreshold_T,screwcrss_coefficient(1)*screwthreshold_T+screwcrss_coefficient(2),'LineWidth',2)

        end
        hold off
        xlabel('Temperature [K]', 'FontSize', 20)
        ylabel('threshold RSS [GPa]', 'FontSize', 20)
    end
    
    % friction coefficient change wrt T
    subplot(2,3,5)
    hold on
    plot(screw_sampling_temperature,screwfriction_coefficient_gpa,'LineWidth',2)
    scatter(screw_sampling_temperature,screwfriction_coefficient_gpa,'r')
    plot(screw_sampling_temperature,mobility_fitting_coefficients(1)*screw_sampling_temperature.^2+mobility_fitting_coefficients(2)*screw_sampling_temperature+mobility_fitting_coefficients(3))
    hold off
    xlabel('Temperature [K]', 'FontSize', 20)
    ylabel('mobility []', 'FontSize', 20)
    
end
    
    Bt_edge_GPa = edgefriction_coefficient_gpa;
    Bt_screw_GPa = mobility_fitting_coefficients(1)*T.^2+mobility_fitting_coefficients(2)*T+mobility_fitting_coefficients(3);
    debye_frequency = debye_frequency_alloy;
    % need T as entry
    edgecrss = edgecrss_coefficient(1)*T+edgecrss_coefficient(2);
    screwcrss = screwcrss_coefficient(1)*T+screwcrss_coefficient(2);
    
% end

