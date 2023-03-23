% loop through each stresses

clear;clc;close all;

input_HEA();

%Dislocation Dynamics simulation in 3-dimension
%
%the best way to run it is by "run dd3d"
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


% links: first 2 collumns-> indices of rn that connects with its neighboring collumn
%          3,4,5 collumns are 'burgers vector'
%          7,8,9 collumns are ?
% rn: first 3 collumns are x,y,z coordinates, 4th collumn is flag

% format longEng
format short

mex -O SegSegForces.c
mex -O SegSegForcesVector.c

if ~exist('restrictSurfaceNodes') %If ~exist, create a bool that allows surface ndoes
    restrictSurfaceNodes=false;
end

if ~exist('writeMovie') %Create the writeMovie variable if not predefined
    writeMovie = false;
end

if writeMovie %Record a movie if specified to do so   
    newVid = VideoWriter('video', 'MPEG-4'); % New
    newVid.FrameRate = 10;
    newVid.Quality = 100;
    open(newVid);
    frame =1;   
end




%% initialization
frame_counter=1;
eventsequence = [];
% eventsequence: 1 straight segment glide
%                2 jog advance in x
%                3 jog expand in z
%                4 create new super jog on straight segment
%                5 connection movement (unsolved)   
%                6 straight segment between two jogs glide


% initialize plots
% f.Name = 'simulation';
% f.Position = [100 100 1400 800];
% set(gcf, 'Position', [10 10 900 800]);
appliedstress_zero = [0 0 0
                        0 0 0
                        0 0 0];

% vn: [A/s]
% unit for time is [s]
% fn fseg: [GPa*A^2]
% Bt: [GPa*s]


for idx_stress = 1:size(stress0,2)
    disp(strcat('stress=',num2str(stress0(idx_stress)),'N/A^2'))
    stress = stress0(idx_stress);

    appliedstress = [0 stress 0
                    stress 0 0
                    0 0 0];   % [N/A^2]
    
    rn = rn_init;
    links = links_init;
    typelist = typelist_init;

    fileID_dtkmc = fopen(strcat('dt_KMC_',num2str(stress0(idx_stress)*1e14),'stress_',num2str(T),'T_',num2str(max_totalDDsteps),'step.txt'),'w');
    fileID_strain = fopen(strcat('strain_time_',num2str(stress0(idx_stress)*1e14),'stress_',num2str(T),'T_',num2str(max_totalDDsteps),'step.txt'),'w');
    fileID_strain_time_kmcjump = fopen(strcat('strain_time_kmcjump',num2str(stress0(idx_stress)*1e14),'stress_',num2str(T),'T_',num2str(max_totalDDsteps),'step.txt'),'w');

    strainP = zeros(max_totalDDsteps,1);
    stress_list = zeros(max_totalDDsteps,1);
    strainP_cumulative = zeros(max_totalDDsteps,1);
    strainP_cumulative = zeros(max_totalDDsteps,1);
    
    DDtime_perstep_list=zeros(max_totalDDsteps,1);
    DDtime_list=zeros(max_totalDDsteps,1);
    totalDDsteps = 0;
    totaltimeline = 0;
    dt_KMC_list = [];

    % cleanup the empty node and link entries at the end of the initial data structures
    [rn,links]=cleanupnodes(rn,links);
   
    % genererate the connectivity list from the list of links
    [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
    consistencycheck(rn,links,connectivity,linksinconnect,rntol);

    % relax the structure first 
    DDtime=0;
    curKMCstep=0;
    for iter_relax = 1:10
        [rn,vn,dt,fn,fseg,DDtime]=feval(integrator,rn,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility,dopartials,stackforcevec,DDtime,rann,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,curKMCstep,printfreq,doshielding,SurfacePlane,dobox,boxD,Bt_edge,Bt_screw,edgecrss,screwcrss); 
    end

%%  KMC loop
for curKMCstep=1:max_totalKMCsteps  % total steps of KMC events followed by DD 
    if totalDDsteps >= max_totalDDsteps
        break
    end
    if totaltimeline >= max_totaltime
        break
    end

    disp(strcat('curstep=',num2str(curKMCstep)))
    if writeMovie && curKMCstep==1 %Open the movie file for writing
        open(newVid);
    end
    
    %   calculate resolved shear stress tau_rss for each segment
    [fseg_for_rss,fpk,fr0,fs0]=segforcevec_debug(MU,NU,a,Ec,rn,links,appliedstress,[],mobility,dopartials,stackforcevec,doinclusion,inclusion_pos_rad,doSFT,SFT_plane,doshielding,SurfacePlane,dobox,boxD);
    %     fseg  in [N]

    %% fieldpoint stress calculation
    idx_jog = [];
    for idx_seg = 1:size(links,1)
        if typelist(idx_seg)==2
            idx_jog = [idx_jog;idx_seg];
        end
    end
    
    if size(idx_jog)>0 % if there is superjog, calculate work first, then calculate rate, finally sample the event
        % generate new typelist because of remesh enabled    
        %%  calculate work and elastic energy change due to migration
        work_list = zeros(size(links,2),4);
        Eelastic_list = zeros(size(links,1),4);% each jog partial has an elastic energy, composed of four terms: forward left and right, backward left and right
    
        for idx_seg = 1:size(links,1)
    %         for each superjog parital, move the entire superjog to calculate
    %         elastic energy change, so the consequence is the total rate is
    %         dobouled, each superjog has two partials (rates)
            if typelist(idx_seg)==2
                idx_r1 = links(idx_seg,1);
                if rn(idx_r1,4) == 7 && rn(idx_r1+1,4) == 7 && rn(idx_r1+2,4) == 7 && rn(idx_r1+3,4) == 7
    %                 left part of super jog
                    idx_jog1 = idx_r1;
                    idx_jog2 = idx_r1+1;
                    idx_jog3 = idx_r1+2;
                    idx_jog4 = idx_r1+3;
                elseif rn(idx_r1,4) == 7 && rn(idx_r1+1,4) == 7 && rn(idx_r1-1,4) == 7 && rn(idx_r1-2,4) == 7
    %                 left part of super jog
                    idx_jog1 = idx_r1-2;
                    idx_jog2 = idx_r1-1;
                    idx_jog3 = idx_r1;
                    idx_jog4 = idx_r1+1;
                else
                    disp('problem with superjog')
                end
    
                %% calculate work
                for iter_direction = 1:2
                    r1= rn(idx_jog1,1:3);
                    r2= rn(idx_jog2,1:3);
                    r3= rn(idx_jog3,1:3);
                    r4= rn(idx_jog4,1:3);
        
                    x = 0.5*(r2+r3);
        
                    seg_start_list = rn(1:end-1,1:3);
                    seg_end_list = rn(2:end,1:3);
                    mu=MU;
                    nu=NU;
                    % Stress(i,:)=[s_11 s_22 s_33 s_12 s_23 s_13]
                    fieldPointStress = FieldPointStress(x,seg_start_list,seg_end_list,b,a,mu,nu); % same as [mu] -> N/A^2
                    tau = [fieldPointStress(1) fieldPointStress(4) fieldPointStress(6);
                            fieldPointStress(4) fieldPointStress(2) fieldPointStress(5);
                            fieldPointStress(6) fieldPointStress(5) fieldPointStress(3)]+appliedstress;
        
                    % calculate RSS and subtrate CRSS
                    l = norm(r3-r2);
                    t = (r3-r2)/l; % line direction unit vector
                    g = cross(n,t);  % glide direction unit vector
                    rss = n*tau*g';
                    stress_available = rss-edgecrss;
                    if stress_available<0
                        stress_available = 0;
                    end
    
                    if iter_direction == 1
                        work_list(idx_seg,1) = -stress_available*w*norm(b)*norm(b);  % forward the barrier is lowered by applied stress
                    else
                        work_list(idx_seg,2) = +stress_available*w*norm(b)*norm(b);
                    end
    
                %% calculate Eelastic 
                    rn_after_migration = rn;
                    if iter_direction == 1
                        rn_after_migration(idx_jog1,:) = rn_after_migration(idx_jog1,:)+[norm(b) 0 0 0];
                        rn_after_migration(idx_jog2,:) = rn_after_migration(idx_jog2,:)+[norm(b) 0 0 0];
                        rn_after_migration(idx_jog3,:) = rn_after_migration(idx_jog3,:)+[norm(b) 0 0 0];
                        rn_after_migration(idx_jog4,:) = rn_after_migration(idx_jog4,:)+[norm(b) 0 0 0];
                    else
                        rn_after_migration(idx_jog1,:) = rn_after_migration(idx_jog1,:)-[norm(b) 0 0 0];
                        rn_after_migration(idx_jog2,:) = rn_after_migration(idx_jog2,:)-[norm(b) 0 0 0];
                        rn_after_migration(idx_jog3,:) = rn_after_migration(idx_jog3,:)-[norm(b) 0 0 0];
                        rn_after_migration(idx_jog4,:) = rn_after_migration(idx_jog4,:)-[norm(b) 0 0 0];
                    end
    
        %           calculating Eelastic for entire superjog
                    fseg_after_migration=segforcevec(MU,NU,a,Ec,rn_after_migration,links,appliedstress,[],mobility,dopartials,stackforcevec,doinclusion,inclusion_pos_rad,doSFT,SFT_plane,doshielding,SurfacePlane,dobox,boxD);
                    
                    fseg_r1side = fseg_after_migration(idx_jog1-1,4:6)+fseg_after_migration(idx_jog1,1:3);
                    z_distance_r1side = rn(idx_jog1,3)-rn(idx_jog1-1,3);
                    Area1 = 0.5*(z_distance_r1side*norm(b));
                    Ee_r1side = abs(fseg_r1side(1)*Area1/norm(b)); % N*A
                    Ee_ev_r1side = Ee_r1side*1e-10*6.242e+18;
        
                    fseg_r2side = fseg_after_migration(idx_jog3,4:6)+fseg_after_migration(idx_jog4,1:3);
                    z_distance_r2side = rn(idx_jog4+1,3)-rn(idx_jog4,3);
                    Area2 = 0.5*(z_distance_r2side*norm(b));
                    Ee_r2side = abs(fseg_r2side(1)*Area2/norm(b)); % N*A
                    Ee_ev_r2side = Ee_r2side*1e-10*6.242e+18;
    
                    if iter_direction == 1
                        Eelastic_list(idx_seg,1) = Ee_ev_r1side; % forward left
                        Eelastic_list(idx_seg,2) = Ee_ev_r2side; % forward right
                    else
                        Eelastic_list(idx_seg,3) = Ee_ev_r1side; % backward left
                        Eelastic_list(idx_seg,4) = Ee_ev_r2side; % backward left
                    end
                end
            end
        end
    
        %  calculate total rate
        N = sum(typelist == 2)/2; % number of superjogs on the line 
%         Efv = Evac_sampling(numbins_vac,weights_vac,energies_vac,cumulative_weights_vac);
        
        Emv_list = [];
        % create Emv list for every jog to have their own Emv (same as their counterpart)
        for iter_superjog = 1:2*N
            if mod(iter_superjog,2)==0
                Emv = Emv_list(iter_superjog-1);
            else
                Emv = Emig_sampling(numbins_mig,weights_mig,energies_mig,cumulative_weights_mig);
            end
            Emv_list = [Emv_list;Emv];
        end

        [R_total,eventlist] = calTotalRate(rn,links,typelist,mu_v,Emv_list,kb,T,work_list,Eelastic_list);
        cumulative_ratelist = calCumulativeRate(eventlist);
        
        q_list = [];
        for i = 1:size(eventlist,1)
            q_list = [q_list;{i,eventlist(i).idx_seg,eventlist(i).segtype,eventlist(i).rate,cumulative_ratelist(i),eventlist(i).action}];
        end
    
    %%  sample the events
    % if R_total ~= 0
        rng(100*sum(clock),'philox');
        eta = rand(1,1);

        dt_KMC = -log(eta)/(R_total/2);
        disp(strcat('dt_KMC=',num2str(dt_KMC),'s'))
        dt_KMC_list = [dt_KMC_list,dt_KMC];

        fprintf(fileID_dtkmc,'%e \n',dt_KMC);

        dt0 = dt_KMC;
    
        sampling_rate = R_total*rand(1,1);
        idx_samplingevent = 1;
        milestone_rate = cumulative_ratelist(idx_samplingevent);
        
        while  sampling_rate>milestone_rate
            idx_samplingevent = idx_samplingevent+1;
            milestone_rate = cumulative_ratelist(idx_samplingevent);
        end
        
        idx_seg_sampled = eventlist(idx_samplingevent).idx_seg;
        type_seg_sampled = eventlist(idx_samplingevent).segtype;
        action_seg_sampled = eventlist(idx_samplingevent).action;
        disp(strcat('idx_seg_sampled=',num2str(idx_seg_sampled)))
        
        rn_beforemigration = rn;  % used to calculate plastic strain due to jog migration


        %% output workspace variables

        output_filename = strcat('output_',num2str(stress0(idx_stress)*1e14),'stress_',num2str(T),'T_',num2str(curKMCstep),'step.mat');
        save(output_filename,'rn','links','plim','doinclusion','inclusion_pos_rad','SurfacePlane','typelist','viewangle','plim_x','plim_y','plim_z','includeSuperjogs','rn_init')

        if curKMCstep ~= 1
        if eventlist(idx_samplingevent).action==1
            jump_distance = norm(b);
        elseif eventlist(idx_samplingevent).action==2
            jump_distance = -norm(b);
        end
    
        % superjog move in x through vacancy jump
        disp(strcat('superjog glide in x =',num2str(jump_distance),'A'))
        %       the line tangent of jog segment is in +y direction    
        if rn(idx_seg_sampled,2)<rn(idx_seg_sampled+1,2)
            rn(idx_seg_sampled,:) = rn(idx_seg_sampled,:)+[jump_distance 0 0 0];   
            rn(idx_seg_sampled+1,:) = rn(idx_seg_sampled+1,:)+[jump_distance 0 0 0];  
        
            rn(idx_seg_sampled+2,:) = rn(idx_seg_sampled+2,:)+[jump_distance 0 0 0];   
            rn(idx_seg_sampled+3,:) = rn(idx_seg_sampled+3,:)+[jump_distance 0 0 0]; 
        end
    
        %       the line tangent of jog segment is in -y direction    
        if rn(idx_seg_sampled,2)>rn(idx_seg_sampled+1,2)
            rn(idx_seg_sampled,:) = rn(idx_seg_sampled,:)+[jump_distance 0 0 0];   
            rn(idx_seg_sampled+1,:) = rn(idx_seg_sampled+1,:)+[jump_distance 0 0 0];  
        
            rn(idx_seg_sampled-1,:) = rn(idx_seg_sampled-1,:)+[jump_distance 0 0 0];   
            rn(idx_seg_sampled-2,:) = rn(idx_seg_sampled-2,:)+[jump_distance 0 0 0]; 
        end

        end
        rn_aftermigration = rn;  % used to calculate plastic strain due to jog migration
        dA_migrate = sum(rn_aftermigration(:,1)-rn_beforemigration(:,1))*w;
        strainP_jogmigrate=dislocationdensity*norm(b)*dA_migrate/2/(totallength);

    end


%% run DD 
% remesh+separation+collision
DDtime=0;  % simulation time  
if size(idx_jog,1) == 0
    dt_KMC = 100; % dummy, used to enter the loop 
end

% run DD for time given by KMC (with superjogs), or default time dt (without superjogs)
flag_firstDD = 1;
while DDtime<dt_KMC && totalDDsteps < max_totalDDsteps && totaltimeline < max_totaltime

    rn_beforeDD = rn;
    DDtime_before = DDtime;
%     DDtime is time spent in each DD single run

    dt0=3e-11; % 1ps = 1e-12s
    dt=3e-10;  % this dt was used in original DD, but when KMC is introduced, dt is replaced with poisson variate
    [rnnew,vn,dt,fn,fseg,DDtime]=feval(integrator,rn,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility,dopartials,stackforcevec,DDtime,rann,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,curKMCstep,printfreq,doshielding,SurfacePlane,dobox,boxD,Bt_edge,Bt_screw,edgecrss,screwcrss); 
    totalDDsteps = totalDDsteps+1;
    rn_afterDD = rnnew;
    DDtime_after = DDtime;
    DDtime_perstep = DDtime_after-DDtime_before;
    disp(strcat('totalDDsteps=',num2str(totalDDsteps)))
    disp(strcat('DDtime_perstep=',num2str(DDtime_perstep),'s'))

    % separation collision and remesh
    rnnew=[rnnew(:,1:3) vn rnnew(:,4)];
    linksnew=links;
    connectivitynew=connectivity;
    linksinconnectnew=linksinconnect;
    fsegnew=fseg;
    stackforcevecnew=stackforcevec;
    if(doseparation)
        %spliting of nodes with 4 or more connections
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=separationnew(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,mobility,MU,NU,a,Ec,rann,appliedstress,dopartials,stackforcevecnew,rann,gamma,rntol,maxconnections,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD);
    end
    
    if(docollision)
    %     collision detection and handling
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=collision(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,appliedstress,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD,dt,Bt_edge,Bt_screw);
    end
    
    if(doremesh)
        %remesh
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=remesh(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,appliedstress,mobility,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD,dt,Bt_edge,Bt_screw);      
    end

    rnnew = sortrows(rnnew,3);
    linksnew = generateLinks(rnnew,b,n);


    rn=[rnnew(:,1:3) rnnew(:,7)];
    vn=rnnew(:,4:6);

    [connectivity,linksinconnect]=genconnectivity(rn,linksnew,maxconnections);
    consistencycheck(rn,linksnew,connectivity,linksinconnect,rntol);


    links=linksnew;

%     linksinconnect=linksinconnectnew;
    fseg=fsegnew;
    stackforcevec=stackforcevecnew;
    typelist = generateNewTypelist(rn,links);
    
    %%%%%%%% only calculate strainP in specific range
    %% store run time information
    DDtime_list=DDtime;  % delta_t is the time defined in KMC
    DDtime_perstep_list=DDtime_perstep; 
    totaltimeline =totaltimeline+DDtime_perstep_list;
    disp(strcat('totaltimeline=',num2str(totaltimeline),'s'))
    disp(' ')

    Area_after = 0;
    Area_before = 0;
    
    % make sure the four vertices for the two superjogs are correctly
    % sorted
    
    rn_beforeDD = sortrn(rn_beforeDD);
    rn_afterDD = sortrn(rn_afterDD);
    
    for iter_seg = 1:size(rn_afterDD,1)-1
        Area_after = Area_after+rn_afterDD(iter_seg,1)*(rn_afterDD(iter_seg+1,3)-rn_afterDD(iter_seg,3)); % first term is distance traveled, second term is spacing between adjacent nodes
    end
    for iter_seg = 1:size(rn_beforeDD,1)-1
        Area_before = Area_before+rn_beforeDD(iter_seg,1)*(rn_beforeDD(iter_seg+1,3)-rn_beforeDD(iter_seg,3)); % first term is distance traveled, second term is spacing between adjacent nodes
    end

    dA = Area_after-Area_before;

    strainP= dislocationdensity*norm(b)*dA/2/(totallength);
    
    strainP_cumulative_direct = dislocationdensity*norm(b)*Area_after/2/(totallength);
    if totalDDsteps==1
        strainP_cumulative=strainP;
    else
        strainP_cumulative=strainP_cumulative+strainP;
    end
    
    strainP_cumulative = strainP_cumulative_direct;

%   first round DD step, print time and plastic strain
    if size(idx_jog,1) ~= 0
        if flag_firstDD == 1
            fprintf(fileID_strain_time_kmcjump,'%8.4f \t',totaltimeline*1e12);
            fprintf(fileID_strain_time_kmcjump,'%8.4e \t  \n',strainP_cumulative);
        end
    end
    flag_firstDD=0;

    stress_list(totalDDsteps)=stress;

    fprintf(fileID_strain,'%8.4f \t',totaltimeline*1e12);
    fprintf(fileID_strain,'%8.4e \t  \n',strainP_cumulative);


    %% write to movie
    plotnodes_alan_normalizedscale(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,typelist,viewangle,plim_x,plim_y,plim_z,includeSuperjogs); 
    drawnow

    if writeMovie && mod(totalDDsteps,20)==1
        theframe = getframe(gcf);
        framelist(frame_counter) = theframe;
        frame_counter=frame_counter+1;
    end
    pause(0.01);
    cla

%     %% output workspace variables
%     steps_to_saveconfig = 1:ceil(max_totalDDsteps/10):max_totalDDsteps;
%     if any(steps_to_saveconfig== totalDDsteps) 
%         output_filename = strcat('output_',num2str(stress0(idx_stress)*1e14),'stress_',num2str(T),'T_',num2str(totalDDsteps),'step.mat');
%         save(output_filename,'rn','links','plim','doinclusion','inclusion_pos_rad','SurfacePlane','typelist','viewangle','plim_x','plim_y','plim_z','includeSuperjogs','rn_init')
%     end

    if size(idx_jog,1) == 0 % for straight dislocation, only one DD run per KMC step
        break; 
    end
     % for superjogged dislocation, run sufficient times of DD, totalling dt_KMC 
        % maybe need to consider one DD step exceeding dt_KMC greatly
    
end

% after DD finished, update DD time into totaltime
disp(' ')

end  % end of KMC loop


%% plot the final configuration
% plotnodes_alan_normalizedscale(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,typelist,viewangle,plim_x,plim_y,plim_z,includeSuperjogs); 
% drawnow

%% output video
if writeMovie %Record a movie if specified to do so   
    newVid = VideoWriter('video', 'MPEG-4'); % New
    newVid.FrameRate = 10;
    newVid.Quality = 100;
    open(newVid);
    for frame = 1:length(framelist)
        writeVideo(newVid,framelist(frame));
    end
    close(newVid);
end


end  % end of looping through stress


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%