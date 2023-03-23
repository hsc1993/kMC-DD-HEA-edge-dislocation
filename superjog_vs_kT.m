clc;clear;close all




numiteratioin = 20;
includeSuperjogs = 1;
flag_plot=0;


doinclusion=0;
Ri=5; % Radius inclusion - to be varied 
inclusion_pos_rad=[0 8 0 Ri 0]; % 11 -11
SurfacePlane=[0 0 0 0];


mean_numvac_list = [];
std_numvac_list = [];



T_low = 300;      % Temperature low end
T_high = 1500;      % Temperature high end

T_step = 100;
num_T = (T_high-T_low)/100+1;


[a0,mu_v] = debye_frequency();


b_globalframe = a0/2*[1 1 1];
b_globalframe_norm = a0/2*norm([1 1 1]);
Va = a0^3/2;   % BCC atomic volume  A^3
w = a0*sqrt(6)/3;  %  1/6[112]

b = [  1 0 0 ]*b_globalframe_norm;
n = [  0 1 0 ];

kb = 8.617333262145e-5;	% eV * K^-1


% generate vacancies along the dislocation line --Alan
segmentlength = sqrt(6)*a0/3;
numdislocationsegments = 1000;
dislocationtotallength = numdislocationsegments*segmentlength;
idx_rn = 1;
yheight = sqrt(2)/2*a0;


f = figure();
f.Name = 'vac and mig histogram';
f.Position = [100 100 1400 800];
subplot(1,2,1)
[numbins_vac,weights_vac,energies_vac,cumulative_weights_vac] = Evac_data_histing();
subplot(1,2,2)
[numbins_mig,weights_mig,energies_mig,cumulative_weights_mig] = Emig_data_histing();


numvac = [];

T_list = [];

for iter_T = 1:num_T
    T = T_low+(iter_T-1)*T_step;
    T_list = [T_list,T]
    disp(strcat('T=',num2str(T)))



for iter = 1:numiteratioin

typelist = [];

idx_vacancy_list = [];
if includeSuperjogs==1
%     dislocation_generation_range = [floor(numdislocationsegments/3),ceil(numdislocationsegments-floor(numdislocationsegments*1/3))];
    dislocation_generation_range = [1,numdislocationsegments];
for i = dislocation_generation_range(1):dislocation_generation_range(2)
    [Efv,Efv_mean] = Evac_sampling(numbins_vac,weights_vac,energies_vac,cumulative_weights_vac);
    pfv = exp(-Efv/(kb*T)); % vacancy formation probability

    p = rand;
    if p<pfv
        idx_vacancy = randi([dislocation_generation_range(1) dislocation_generation_range(2)]);
        if any(idx_vacancy_list==idx_vacancy)
            % do not generate vacancy when there is one already
            continue
        end
%         if any(idx_vacancy_list==idx_vacancy-1) || any(idx_vacancy_list==idx_vacancy+1)
%             % do not generate vacancy near an existing vacancy
%             continue
%         end
        idx_vacancy_list = [idx_vacancy_list,idx_vacancy];
    end
end
idx_vacancy_list = sort(idx_vacancy_list);


end
flag_vac = 0;

rn = [0 0 floor(1-numdislocationsegments/2)*segmentlength 0];   % first node position

dislocation = dislocation;

for i = 1:numdislocationsegments
    z_multiplyer = floor(i-numdislocationsegments/2);

    if flag_vac == 1
        flag_vac = 0;
        continue
    end
    if any(idx_vacancy_list==i)
        rn_jog1 = [0 0 (z_multiplyer)*segmentlength 7];
        rn_jog2 = [0 yheight (z_multiplyer)*segmentlength 7];
        
        rn_jog3 = [0 yheight (z_multiplyer+1)*segmentlength 7];
        rn_jog4 = [0 0 (z_multiplyer+1)*segmentlength 7];
        
        rn_newstraight = [0 0 (z_multiplyer+2)*segmentlength 0];
        rn = [rn(1:end-1,:);rn_jog1;rn_jog2;rn_jog3;rn_jog4;rn_newstraight];
        
        idx_rn = idx_rn+4;
        typelist = [typelist(1:end-1);3;2;1;2;3];
        flag_vac = 1;


%         create jog objects
        jog1 = segment;
        jog2 = segment;
        jog1.idx_pairing_jog = i+1;
        jog2.idx_pairing_jog = i;
        jog1.idx_seg = i;
        jog2.idx_seg = i+1;
        jog1.type = 'jog';
        jog2.type = 'jog';
        jog1.rn_start = [0 0 (z_multiplyer)*segmentlength 7];
        jog1.rn_end = [0 yheight (z_multiplyer)*segmentlength 7];
        jog2.rn_start = [0 yheight (z_multiplyer+1)*segmentlength 7];
        jog2.rn_end = [0 0 (z_multiplyer+1)*segmentlength 7];
        
        newstraight = segment;
        newstraight.idx_seg = i+2;
        newstraight.type = 'straight';
        newstraight.idx_pairing_jog = 'None';
        newstraight.rn_start = [0 0 (z_multiplyer+2)*segmentlength 0];
        newstraight.rn_end = [0 0 (z_multiplyer+3)*segmentlength 0];

        superjog = superjog;
        superjog.idx_seg_list = [[0 0 (z_multiplyer)*segmentlength 7];
            [0 yheight (z_multiplyer)*segmentlength 7];
            [0 yheight (z_multiplyer+1)*segmentlength 7];
            [0 0 (z_multiplyer+1)*segmentlength 7]];

        dislocation.segmentlist = dislocation.segmentlist(1:end-1);
%         dislocation.segmentlist = [dislocation.segmentlist;jog1;jog2;newstraight];
        dislocation.segmentlist{end+1} = superjog;

    else
        rn_straight = [0 0 (z_multiplyer+1)*segmentlength 0];
        rn = [rn;rn_straight];
%         links = [links;[idx_rn idx_rn+1 b n]];
        idx_rn = idx_rn+1;
        typelist = [typelist;1];

%         create straight segment objects
        straightseg = segment;
        straightseg.idx_seg = i;
        straightseg.type = 'straight';
        straightseg.idx_pairing_jog = 'None';
        straightseg.rn_start = [0 0 (z_multiplyer)*segmentlength 0];
        straightseg.rn_end = [0 0 (z_multiplyer+1)*segmentlength 0];

%         dislocation.segmentlist = {dislocation.segmentlist;straightseg};
        dislocation.segmentlist{end+1} = straightseg;
    end
end

rn(1,4) = 7;
rn(end,4) = 7;
rn_init = rn;


links = generateLinks(rn,b,n);




plim_Y = 2*max(abs(rn(:,2)));
if plim_Y == 0
    plim_Y = yheight*2;
end
plim_X = 10*plim_Y*1;
plim_Z_end = rn(end,3)*1;
plim_Z_1 = rn(1,3)*1;

plim_x = [plim_Z_end,plim_Z_1];
plim_y= plim_X;
plim_z= plim_Y;

plim = 1;  % dummy, not used in plotnodes_alan

rmax=0.1;
 viewangle = [1,1,1]; %Upward tilted - able to see top and climb


if flag_plot==1
    f = figure();
    f.Name = 'Initial Configuration';
    f.Position = [100 100 1000 800];
    plotnodes_alan_normalizedscale(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,typelist,viewangle,plim_x,plim_y,plim_z,includeSuperjogs); 
end

numvac = [numvac;size(idx_vacancy_list,2)/dislocationtotallength*10];
mean_numvac = mean(numvac);
std_numvac = std(numvac);

end


mean_numvac_list = [mean_numvac_list,mean_numvac];
std_numvac_list = [std_numvac_list,std_numvac];

end

close all

kt_inverse = 1/kb./T_list;
hold on
plot(kt_inverse,mean_numvac_list,'LineWidth',4)
% scatter(kt_inverse,mean_numvac_list,20)
N0 = floor(dislocationtotallength/w)/dislocationtotallength*exp(-(Efv_mean.*kt_inverse))*10;   % *10 to calculate density per nm

plot(kt_inverse,N0,'LineWidth',4)

lg = legend('actual number of vacancies','expected number of vacancies using avarage E^f_v');
lg.FontSize = 14;
hold off

title(strcat('number of vacancies on a ', num2str(numdislocationsegments) ,'b long straight edge dislocation'),'FontSize', 30)
xlabel('1/kT','FontSize', 30)
ylabel('N','FontSize', 30)
set(gca, 'YScale', 'log')




% mean vacancy 
namestring = 'mean_vacancy.txt';
Fid = fopen(namestring,'w');

for j = 1:size(kt_inverse,2)
    fprintf(Fid,'%i  %i',kt_inverse(j),mean_numvac_list(j));
    fprintf(Fid,'\n');
end

% mean vacancy 
namestring = 'std_vacancy.txt';
Fid = fopen(namestring,'w');

for j = 1:size(kt_inverse,2)
    fprintf(Fid,'%i  %i',kt_inverse(j),std_numvac_list(j));
    fprintf(Fid,'\n');
end


% expected mean vacancy 
namestring = 'N0_vacancy.txt';
Fid = fopen(namestring,'w');

for j = 1:size(kt_inverse,2)
    fprintf(Fid,'%i  %i',kt_inverse(j),N0(j));
    fprintf(Fid,'\n');
end









