%Dislocation Dynamics simulation in 3-dimension
%
%the best way to run it is by "rundd3d"
%
%Features:
%infinite boundary condition (no pbc)
%linear mobility law (mobfcc0,mobfcc1)
%no remesh, no collision detection, no reaction
%N^2 interaction (no neighbor list, no fast multipole)

%Data structure:
%NMAX:    maximum number of nodes (including disabled ones)
%LINKMAX: maximum number of links (including disabled ones)
%rn: (NMAX,4) array of nodal positions (last column is flag: -1 means disabled)
%vn: (NMAX,3) array of nodal velocities
%fn: (NMAX,3) array of nodal forces
%links: (LINKMAX,8) array of links (id1,id2,bx,by,bz,nx,ny,nz)

%default value if run by itself (e.g. not through "rundd3d")
if(~exist('rn'))
    initparams;  %default parameter settings (can be override by restart files below)
    load frsource; %read in dislocation configuration and setting files
end
if(~exist('dt'))
    dt=dt0;
end
if(~exist('doshielding'))
    doshielding=0;
end
if(~exist('doSFT'))
    doSFT=0;
end
if(~exist('SFT_plane'))
    SFT_plane=[];
else
    SFTtol=1e-5;
    for i=0:doSFT-1
        SFT_index=4*i+1;
        SFT_plane(SFT_index,4:6)=SFT_plane(SFT_index,4:6)-SFTtol.*SFT_plane(SFT_index,1:3);
        SFT_plane(SFT_index+1,4:6)=SFT_plane(SFT_index+1,4:6)-SFTtol.*SFT_plane(SFT_index+1,1:3);
        SFT_plane(SFT_index+2,4:6)=SFT_plane(SFT_index+2,4:6)-SFTtol.*SFT_plane(SFT_index+2,1:3);
        SFT_plane(SFT_index+3,4:6)=SFT_plane(SFT_index+3,4:6)-SFTtol.*SFT_plane(SFT_index+3,1:3);
    end
end
if(~exist('doinclusion'))
    doinclusion=0;
    inclusion_pos_rad=0;
end

% cleanup the empty node and link entries at the end of the initial data structures
[rn,links]=cleanupnodes(rn,links);
%plot dislocation structure
figure(1); 
plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad); view(viewangle); xlim([-plim plim]); ylim([-plim plim]); zlim([-plim plim]);
drawnow

%add the partials to the code. A flag is needed to have the right sense of the stacking fault force(Enrique Sep 2005)


if((findstr(mobility,'fcc'))&(dopartials))
    numlinks=size(links,1);
    if(~exist('stackforcevec'))
        stackforcevec=zeros(numlinks,3);
    end
    %[connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
    %consistencycheck(rn,links,connectivity,linksinconnect,rntol);
    %[rn,links,stackforcevec]=stackfaultsplit(rn,links,connectivity,rann,gamma,stackforcevec,maxconnections);
end

% cleanup the empty node and link entries at the end of the initial data structures
[rn,links]=cleanupnodes(rn,links);



% genererate the connectivity list from the list of links
[connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
consistencycheck(rn,links,connectivity,linksinconnect,rntol);

%plot dislocation structure
figure(1); 
plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad); view(viewangle); xlim([-plim plim]); ylim([-plim plim]); zlim([-plim plim]);
drawnow

data=zeros(totalsteps,1);

totaltime=0;
dt=min(dt,dt0);
mdold=10;
counter=0;

%%%%% Parameters for Stress Strain Curve

glideplane=[1 1 1]/norm()
for i=1:num_tot_points
    
for curstep=1:totalsteps,
    
    %integrating equation of motion
    [rnnew,vn,dt,fn,fseg,totaltime]=feval(integrator,rn,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility,dopartials,stackforcevec,totaltime,rann,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,curstep,printfreq,doshielding);
    % plastic strain and plastic spin calculations
    [ep_inc,wp_inc]=calcplasticstrainincrement(rnnew,rn,links,(2*plim)^3);
    
    numprint=length(printnode);
    vnmax=size(vn,1);
    
    if(mod(curstep,printfreq)==0)
        disp(sprintf('step%3d dt=%e\ttime=%e',...
            curstep,dt,totaltime));
        for h=1:numprint
            if(printnode(h)>vnmax)
                printnode(h)=vnmax-(h-1);
            end
            disp(sprintf('v%d=(%e,%e,%e)\nforce%d=(%e, %e, %e)\n',...
            printnode(h),vn(printnode(h),1),vn(printnode(h),2),vn(printnode(h),3),printnode(h),fn(printnode(h),1),fn(printnode(h),2),fn(printnode(h),3)));
        end
    end
    if(mod(curstep,plotfreq)==0)
        %figure(1); 
        %plim=max(max(abs(rn(:,1:3))))+20;
        plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad);  xlim([-plim plim]); ylim([-plim plim]); zlim([-plim plim]);
        view(viewangle);
        drawnow
        pause(0.01);
    end
    
    rnnew=[rnnew(:,1:3) vn rnnew(:,4)];
    linksnew=links;
    connectivitynew=connectivity;
    linksinconnectnew=linksinconnect;
    fsegnew=fseg;
    stackforcevecnew=stackforcevec;
    if(docrossslip)
        %Fleischer and Escaig Cross-slip
        crossslip=0;
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew,crossslip,printnode]=crossslipfuncnewFrank(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,crossslip,integrator,dt,dt0,rmax,totaltime,areamin,lmin,printnode,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding);
    end
    if(doseparation)
        %spliting of nodes with 4 or more connections
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=separationnew(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,mobility,MU,NU,a,Ec,rann,appliedstress,dopartials,stackforcevecnew,rann,gamma,rntol,maxconnections,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding);
    
    end
    if(docollision)
        %collision detection and handling
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=collision(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,appliedstress,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding);
    end
    
    if(doremesh)
        %remesh
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=remesh(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,appliedstress,mobility,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding);      
    end
   
    rn=[rnnew(:,1:3) rnnew(:,7)];
    vn=rnnew(:,4:6);
    links=linksnew;
    connectivity=connectivitynew;
    linksinconnect=linksinconnectnew;
    fseg=fsegnew;
    stackforcevec=stackforcevecnew;
    %store run time information
    %time step
    data(curstep,1)=dt;
    save restart
    consistencycheck(rn,links,connectivity,linksinconnect,rntol);
    conservstakforce(stackforcevec,links,rn,connectivity);
   
end
ep_inc
wp_inc
save restart

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%