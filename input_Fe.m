%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input file for dislocation - precipitate interaction   %  
% Eshelby inclustions                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% From Ninive et al: 
% - pure Al GGA-PBE
% mu = 4.5628e-09 N/A2*a2, E = 6.38793e-09 N/A2*a2, B = -1.16291e+08 (N/A2*a2)^-1

writeMovie = true; %Bool to save a movie or not
restrictSurfaceNodes = false; %Bool to keep nodes off of the particles surface or not


sim_params = simulation_parameters; %Define all the simulation parameters needed for emission/vacancies

% Matrial constants Al
gamma=6.612e-11; %Lu et al 1999
MU = 4.5628e-09; % Ninive 2014
%MU = 5.28e-9; %Liu 2004
NU = 0.332;% Ninive 2014
a=1/sqrt(2);%b %1.75*(1/sqrt(2));% 3b  lmin/sqrt(3)*0.5;
%Ec = (1/(pi*(a/5)^2)).*[1.897e-10 1.135e-10 6.91e-11 5.195e-10 6.57e-10 6.57e-10];% [Perfect Shockley StairRod Hirth Frank] [J/l3]MU/(4*pi)*log(a/0.1);
%Ec = [MU/(4*pi)*log(a/0.1) (1/(pi*(a/5)^2))*1.135e-10 (1/(pi*(a/5)^2))*6.91e-11 (1/(pi*(a/5)^2))*5.195e-10 (1/(pi*(a/5)^2))*6.57e-10 (1/(pi*(a/5)^2))*6.57e-10];% [Perfect Shockley StairRod Hirth Frank] [J/l3]MU/(4*pi)*log(a/0.1);
Ec = [MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1) MU/(4*pi)*log(a/0.1)];% [Perfect Shockley StairRod Hirth Frank] [J/l3]MU/(4*pi)*log(a/0.1);

Ri=5; % Radius inclusion - to be varied 

% % Initial dislocation structure edge
rn    = [
         0 0 -60   7;
         0 0 0  0;
         0 0 60   7;        
     ];

b1 = [  0 1 0 ]/sqrt(2);
n1 = [  1 0 0 ];


links = [
         1   2  b1 n1;
         2   3  b1 n1;
       ];

sigma = -50.257e-11.* [ 0  1  0
         1  0 0
         0 0 0 ];
     
appliedstress = sigma;   

% Inclusion
doinclusion=1;


inclusion_pos_rad=[0 8 0 Ri 0]; % 11 -11
                  %-200  200 0 30 2.66e-9 0.3045 1e-9 Biscrew Biedge;
                   %xc yc zc R c_1]
                   
%coordinate system 1 (cubic coordinate system)
% e1p = [1 0 0]; e2p = [0 1 0]; e3p = [0 0 1];
% e1p=e1p/norm(e1p); e2p=e2p/norm(e2p); e3p=e3p/norm(e3p);            
% %coordinate system 2
% %rotate so:
% % - x along direction dislocation will move (burgers vector if egde
% % dislocation)
% % - y along glide plane normal
% % - z along dislocation line
% e2  = n1;
% e3 = rn(3,1:3)-rn(1,1:3);
% e1 = cross(e2,e3);
% e1=e1/norm(e1);
% e2=e2/norm(e2); 
% e3=e3/norm(e3);     

maxconnections=50;
lmax = 4;
lmin = 1.25;

rann = a;
rntol = rann/2;
areaminmag(1)=2*lmin*rann*cos(asin(2*rann/(2*pi*lmin))); 
areaminmag(2)=2*lmax*rntol*cos(asin(2*rntol/(2*pi*lmax)));
areamin=min(areaminmag);
areamax=(1/4)*lmax*lmax;

totalsteps=20;
dt0=3e-11;
dt=3e-10;
%mobility='mobfcc0inga';
mobility='mobbcc_climb';
integrator='int_eulerbackward';
%dopartials=0;
doSFT=0;

docrossslip=0;
doremesh=1;
docollision=0;
doseparation=0;
SurfacePlane=[0 0 0 0];
        
%rotation matrix
% Q = [ dot(e1,e1p) dot(e2,e1p) dot(e3,e1p)
%       dot(e1,e2p) dot(e2,e2p) dot(e3,e2p)
%       dot(e1,e3p) dot(e2,e3p) dot(e3,e3p) ];
%   
% %Transform sigma into coordinate system 1
% appliedstress = Q*sigma*Q';
%16.257e-14. = conversion factor for Aluminium
%appliedstress = -1500*13.0676e-14.*(1/(2*sqrt(6))).*([-2 0 -1; 0 2 1; -1 1 0]+(1/3).*[2 0 -1; 0 -2 1; -1 1 0]);



%viewangle=[0,90];  %Looking down z-axis - top view
%viewangle = [-5,10,10]; %From behind the dislocation
%viewangle = [0,10,0]; %Looking straight at the percipitate - climb angle
viewangle = [-8,0,0]; %profile view
%viewangle = [1,3,1]; %Upward tilted - able to see top and climb
%viewangle = [-1,3,1]; %Upward tilted - able to see top and climb


printfreq=10;      
printnode=10;
plotfreq=1;       
plim=30; 
rmax=0.1;
