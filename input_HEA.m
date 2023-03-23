%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input file for dislocation - precipitate interaction   %  
% Eshelby inclustions                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all;

% From Ninive et al: 
% - pure Al GGA-PBE
% mu = 4.5628e-09 N/A2*a2, E = 6.38793e-09 N/A2*a2, B = -1.16291e+08 (N/A2*a2)^-1

writeMovie = true; %Bool to save a movie or not
restrictSurfaceNodes = false; %Bool to keep nodes off of the particles surface or not
sim_params = simulation_parameters; %Define all the simulation parameters needed for emission/vacancies

includeSuperjogs = 1;


% run with jogs 
% totalsteps=5;
% dt0=3e-11; % 1ps = 1e-12s
% dt=3e-10;  % this dt was used in original DD, but when KMC is introduced, dt is replaced with poisson variate


% run pure straight dislocation
max_totalDDsteps=100;
max_totalKMCsteps = 100; % always smaller than or equal max_totalDDsteps
max_totaltime = 30e-12;  % simulation time can not exceed % ps

max_totalDDsteps=5000000;
max_totalKMCsteps = 50000000; % always smaller than or equal max_totalDDsteps
max_totaltime = 30000e-12;  % simulation time can not exceed % ps



dt0=3e-11; % 1ps = 1e-12s
dt=3e-10;  % this dt was used in original DD, but when KMC is introduced, dt is replaced with poisson variate
stress0 = 1 *1e-11; %in units of N/A^2, before exponential is in units of GPa


dorelaxation = 0;
totalsteps_relax = 20;

strainrate=1e-03; % s^-1

% number of segments for dislocation
numSegs = 75;


%% Temperature
T_input = 14;  % 14 means 1400K, 8 means 800K
iter_T = T_input-2;   % 1 to 13
T_low = 300;        % Temperature low end
T_high = 1500;      % Temperature high end
T_step = 100;
num_T = (T_high-T_low)/100+1;

T = T_low+(iter_T-1)*T_step;
disp(strcat('T=',num2str(T)))


%% Matrial constants Al
gamma=6.612e-11; %Lu et al 1999    related to SF
MU = 4.5628e-09; % Ninive 2014    
MU = 94*1e-11; % shear modulus   94 GPa=94e-11 N/A^2
%MU = 5.28e-9; %Liu 2004

NU = 0.33;% Ninive 2014     poisson ratio
a=1/sqrt(2);%b %1.75*(1/sqrt(2));% 3b  lmin/sqrt(3)*0.5;

Ec = [MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1)];% [Perfect Shockley StairRod Hirth Frank] [J/l3]MU/(4*pi)*log(a/0.1);
Ri=5; % Radius inclusion - to be varied 

% - x along direction dislocation will move (burgers vector if egde
% dislocation)
% - y along glide plane normal
% - z along dislocation line

%% Bt and crss extraction
flag_plot = 1;

[a0,mu_v] = debye_frequency();
[Bt_edge, edgecrss] = edge_Bt_crss(T,a0,flag_plot); % edgecrss [N/A^2]  Bt_edge [N/A^2*s]
[Bt_screw, screwcrss] = screw_Bt_crss(T,a0,flag_plot);

Bt_climb = Bt_edge*100;  % make Bt for climb very big to supprese climb from happening

mu_0_prime = mu_v; % vacancy nucleation frequency
% the average Debye frequency of the alloy from the four Debye frequencies of the four constituent elements


%% simulation constants

b_globalframe = a0/2*[1 1 1];
b_globalframe_norm = a0/2*norm([1 1 1]);
Va = a0^3/2;   % BCC atomic volume  A^3
w = a0*sqrt(6)/3;  %  1/6[112]
omega_tensor_v = [-0.02*Va 0 0;
                0 -0.02*Va 0;
                0 0 -0.02*Va];

omega_tensor_m = [0.02 0 0;
                0 0.02 0;
                0 0 0.02];
b = [  1 0 0 ]*b_globalframe_norm;
n = [  0 1 0 ];

kb = 8.617333262145e-5;	% eV * K^-1


%% dislocation line configuration
segmentlength = sqrt(6)*a0/3;
totallength = numSegs*segmentlength;
fraction_dislocation_w_superjogs = 3;
dislocationdensity = 5e-8; % A^-2


%% formation energy and migration energy
% f.Name = 'Evac and Emig histogram';
% f.Position = [100 100 1400 800];
[numbins_vac,weights_vac,energies_vac,cumulative_weights_vac] = Evac_data_histing();
% title('Evac')
[numbins_mig,weights_mig,energies_mig,cumulative_weights_mig] = Emig_data_histing();
% title('Emig')

[~,Efv_mean] = Evac_sampling(numbins_vac,weights_vac,energies_vac,cumulative_weights_vac);
N0 = floor(totallength/w)*exp(-Efv_mean/(kb*T));


maxnumjogs = 1;
%% vacancy formation sampling
% if pfv is used, that is to sample vacancy formation, else it is to insert
% vacancy by brute force at specific location
idx_vacancy_list = [];
if includeSuperjogs==1
    
    for i = 3:numSegs-3
        [Efv,Efv_mean] = Evac_sampling(numbins_vac,weights_vac,energies_vac,cumulative_weights_vac);
        pfv = exp(-(Efv)/(kb*T)); % vacancy formation probability
        p = rand;
%         if p<pfv 
%           if 1<0
%         if i== ceil(numSegs/6) || i== ceil(numSegs/2) || ceil(numSegs*5/6)   % 1/2
%         if i== ceil(numSegs/6) || i== ceil(numSegs/2) || i==ceil(5*numSegs/6)   % 1/2
        if i== ceil(numSegs/3) || i== ceil(2*numSegs/3)   % 1/3 etc
%         if i== ceil(numSegs/4) || i== ceil(2*numSegs/4) || i== ceil(3*numSegs/4)  % 1/4 etc 
%         if i== ceil(numSegs/5) || i== ceil(2*numSegs/5) || i== ceil(3*numSegs/5) || i== ceil(4*numSegs/5)   % 1/5 etc
%         if i== ceil(numSegs/6) || i== ceil(2*numSegs/6) || i== ceil(3*numSegs/6) || i== ceil(4*numSegs/6) || i== ceil(5*numSegs/6)   % 1/6 etc
            idx_vacancy = i;
            idx_vacancy_list = [idx_vacancy_list,idx_vacancy];
        end
    end
    idx_vacancy_list = sort(idx_vacancy_list);
end


rn = [0 0 floor(1-numSegs/2)*segmentlength 0];   % first node position
yheight = sqrt(2)/2*a0;
idx_rn = 1;
flag_vac = 0;


for i = 1:numSegs
    z_multiplyer = floor(i-numSegs/2);

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
        flag_vac = 1;

    else
        rn_straight = [0 0 (z_multiplyer+1)*segmentlength 0];
        rn = [rn;rn_straight];
        idx_rn = idx_rn+1;

    end
end
rn = rn(3:end-2,:);
rn(1,4) = 7;
rn(end,4) = 7;
rn_init = rn;

links = generateLinks(rn,b,n);
links_init = links;
typelist = generateNewTypelist(rn,links);
typelist_init = typelist;


%    sigma is in N/A^2 
multiplier_GPaToeV = 1/160.2176621;


% Inclusion
doinclusion=0;
inclusion_pos_rad=[0 8 0 Ri 0]; % 11 -11

maxconnections=100;
lmax = 10;
lmin = 1.25;
rann = a;
rntol = rann/2;
areaminmag(1)=2*lmin*rann*cos(asin(2*rann/(2*pi*lmin))); 
areaminmag(2)=2*lmax*rntol*cos(asin(2*rntol/(2*pi*lmax)));
areamin=min(areaminmag);
areamax=(1/4)*lmax*lmax;

%% simulation parameters
% mobility='mobbcc_climb';
mobility='mobbcc_climb_with_Bt_varying';
integrator='int_eulerbackward';

dopartials=0;
stackforcevec=zeros(size(rn,1),3);

doSFT=0;
SFT_plane=[];
doshielding=0;
dobox=0;
boxD = [-202.2 202.2; -202.2 202.2; -202.2 202.2];
docrossslip=0;
doremesh=1;
docollision=1;
doseparation=0;
SurfacePlane=[0 0 0 0];
     
% %% use this part to calculate mechanical work in vacancy formation
% %% and need to change vacancy formation sampling to assist the calculation
% sigma_Efv = -50.257e-11.* [ 0  1  0
%          1  0 0
%          0 0 0 ];
% appliedstress_Efv = sigma_Efv;
% Efv_mechanicalwork = Efv_work(MU,NU,rn,b,a,appliedstress_Efv,omega_tensor_v);%eV

%rotation matrix
% Q = [ dot(e1,e1p) dot(e2,e1p) dot(e3,e1p)
%       dot(e1,e2p) dot(e2,e2p) dot(e3,e2p)
%       dot(e1,e3p) dot(e2,e3p) dot(e3,e3p) ];
%   
% %Transform sigma into coordinate system 1
% appliedstress = Q*sigma*Q';
%16.257e-14. = conversion factor for Aluminium
%appliedstress = -1500*13.0676e-14.*(1/(2*sqrt(6))).*([-2 0 -1; 0 2 1; -1 1 0]+(1/3).*[2 0 -1; 0 -2 1; -1 1 0]);

% viewangle = [1,3,1]; %Upward tilted - able to see top and climb
% viewangle=[-1,1,-0.5];  %Looking down z-axis - top view

%viewangle=[0,90];  %Looking down z-axis - top view
%viewangle = [-5,10,10]; %From behind the dislocation
%viewangle = [0,10,0]; %Looking straight at the percipitate - climb angle
% viewangle = [-8,0,0]; %profile view
% viewangle = [1,3,1]; %Upward tilted - able to see top and climb
 viewangle = [1,1,1]; %Upward tilted - able to see top and climb


% x: glide direction y: superjog height z: dislocation length 
% XYZ: for plot purpose, x y z -> Z X Y
printfreq=10;      
printnode=10;
plotfreq=1;       

 
plim_Y = 2*max(abs(rn(:,2)));
if plim_Y == 0
    plim_Y = yheight*2;
end
plim_X = 3*plim_Y;
plim_Z_end = rn(end,3);
plim_Z_1 = rn(1,3);

plim_x = [plim_Z_end*1.2,plim_Z_1*1.2];
plim_y= plim_X;
plim_z= plim_Y;

plim_x = [-100 100];
plim_y= [-50 50];
plim_z= 10;


plim = 1;  % dummy, not used in plotnodes_alan

rmax=0.1;

f.Name = 'Initial Configuration';
f.Position = [100 100 1000 800];
plotnodes_alan_normalizedscale(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,typelist,viewangle,plim_x,plim_y,plim_z,includeSuperjogs); 




