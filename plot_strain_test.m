clear;clc;close all;

stress = [200 400 600 800 1000 1200 1400 1600 1800 2000];
stress = [ 10000 ];
T = 500;
totalsteps = 10;

sca_list =[];
legend_list = {};
P1_list = [];
P2_list = [];
stress_list = [];
temperature_list = [];

% f = figure();
% f.Name = 'plastic strain vs time';
% f.Position = [100 100 1400 800];


idx_start = [1 1 1 1 1 1 1 1 1 1 ];
idx_end = [0 0 0 0 0 0 0 0 0 0];
% checkidx = 10;
% stress = stress(checkidx);
% idx_start = idx_start(checkidx);
% idx_end = idx_end(checkidx);

for idx_stress = 1:size(stress,2)
fileID = fopen(strcat('test_strain_time_',num2str(stress(idx_stress)),'stress_',num2str(T),'T_',num2str(totalsteps),'step.txt'),'r');
formatSpec = '%f';
sizeA = [2 Inf];

A = fscanf(fileID,formatSpec,sizeA);
B = [0;0];
for i = 1:size(A,2)
    if A(2,i)==0
        continue
    end
    B = [B,[A(1,i);A(2,i)]];
end

totalTime = B(1,:);
strainP_cumulative = B(2,:);
hold on
sca(idx_stress) = plot(totalTime,strainP_cumulative,'LineWidth',5,'DisplayName',strcat('stress=',num2str(stress(idx_stress)),'MPa'));
sca_list = [sca_list,sca(idx_stress)];
legend_list{idx_stress} = strcat('stress=',num2str(stress(idx_stress)),'MPa');


P = polyfit(totalTime(idx_start(idx_stress):end-idx_end(idx_stress)),strainP_cumulative(idx_start(idx_stress):end-idx_end(idx_stress)),1);
scatter(totalTime(idx_start(idx_stress)),strainP_cumulative(idx_start(idx_stress)),'LineWidth',52)
scatter(totalTime(end-idx_end(idx_stress)),strainP_cumulative(end-idx_end(idx_stress)),'LineWidth',52)

yfit = P(1)*totalTime+P(2);
plot(totalTime,yfit,'r-.','DisplayName','');
hold off

P1_list = [P1_list;P(1)];
P2_list = [P1_list;P(2)];
stress_list = [stress_list;stress(idx_stress)];
temperature_list = [temperature_list;T];

end

lg = legend (sca_list,legend_list);
lg.FontSize = 18;

xlabel('time [ps]','FontSize',20)
ylabel('plastic strain ','FontSize',20)
set(gca,'FontSize',20)



%% write lists to txt 
fileID = fopen(strcat('T',num2str(T),'.txt'),'w');
for i = 1:size(stress_list,1)
    fprintf(fileID,'%d %d %e \n',temperature_list(i),stress_list(i),P1_list(i));
end


















