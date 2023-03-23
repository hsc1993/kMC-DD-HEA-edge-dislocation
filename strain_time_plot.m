clear;clc;close all;

stress0 = 2;
T = 1500;
totalsteps = 200;

fileID = fopen(strcat('strain_time_',num2str(stress0),'stress_',num2str(T),'T_',num2str(totalsteps),'step.txt'),'r');
formatSpec = '%f';
sizeA = [2 Inf];

A = fscanf(fileID,formatSpec,sizeA);


f = figure();
f.Name = 'plastic strain vs time';
f.Position = [100 100 1400 800];

totalTime = A(1,:);
strainP_cumulative = A(2,:);
hold on
plot(totalTime*1e12,strainP_cumulative,'LineWidth',5)
scatter(totalTime*1e12,strainP_cumulative,'LineWidth',5)
hold off
xlabel('time [ps]','FontSize',20)
ylabel('plastic strain ','FontSize',20)
