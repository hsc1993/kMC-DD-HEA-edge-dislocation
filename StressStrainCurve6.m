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

format longEng

mex -O SegSegForces.c
mex -O SegSegForcesVector.c

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
if(~exist('dopartials'))
    dopartials=0;
    stackforcevec=[0  0  0;
                   0  0  0;
                   0  0  0;
                   0  0  0];
end
if(~exist('dobox'))
    dobox=0;
    boxD = [-202.2 202.2; -202.2 202.2; -202.2 202.2];
end
%lmax_step=[8,6,4,2];
%step_step=[10,10,10,10];
%dist=[8,16,24,32];
tot_shear=[1e-4,1e-4,2e-4,1e-4,1e-4,1e-4,1e-4];
%le=[25,28,31,34,37,40,43];
le=[50 55 60 65 70 75 80 90 100 110];
cyclesrun=1;
for radius=3:2:13
    for i=1:length(le)
        if i>1 || (i==1 && cyclesrun>1); 
           clear vn rn dt dt0 lmin lmax a rann rntol areaminmag areamin aremax mini connectivity links;
           run Input/input_inclusion_Al_rot.m;
        end

        % % Initial dislocation structure edge
        rn    = [
             0 0 -le(i)   7;
             0 0 0  0;
             0 0 le(i)   7;        
         ];

        %lmin=lmin_step(i);
        %lmax=lmax_step(i);

        %radius=5;
        eps=8e-5;
        total_shear=1e-4;
        stress_step=5;

        %load restart.mat
        if radius ==1 
            continue;
        end
        if radius == 0
            doinclusion=0;
            inclusion_pos_rad=0;
            inclusion_pos_rad(1)=0;
            inclusion_pos_rad(2)=15;
            inclusion_pos_rad(3)=0.0001;
            inclusion_pos_rad(4)=radius;
            inclusion_pos_rad(6)=(radius+0.5);
            appt=0.01166066;
            bppt=0.14993865;
            cppt=appt*radius+bppt; %0.2+radius*0.01; 
            inclusion_pos_rad(5)=cppt;
        else
            doinclusion=1;
            inclusion_pos_rad(1)=0;
            inclusion_pos_rad(2)=15;
            inclusion_pos_rad(3)=0.0001;
            inclusion_pos_rad(4)=radius;
            inclusion_pos_rad(6)=(radius+0.5);
            appt=0.01166066;
            bppt=0.14993865;
            cppt=appt*radius+bppt; %0.2+radius*0.01; 
            inclusion_pos_rad(5)=cppt;
        end
        % Remember to change names!
        filename = sprintf('Stress_Strain_curves/Stress_Strain_curve_r%02.0f_length%03.0f_100_ppt_edge.txt',radius,le(i));
        videoname = sprintf('Stress_Strain_curves/dislocation_r%02.0f_length%03.0f_100_ppt_edge.avi',radius,le(i));
        if exist(filename,'file')
            disp(sprintf('File exists. Simulation is already run.'))
            continue;
        end
        cyclesrun=cyclesrun+1;
        p_file=fopen(filename,'w');
        vidObj = VideoWriter(videoname);
        vidObj.Quality = 100;
        vidObj.FrameRate = 10;
        open(vidObj);

        % cleanup the empty node and link entries at the end of the initial data structures
        [rn,links]=cleanupnodes(rn,links);
        %plot dislocation structure
        %figure(1); 
        %plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad); view(viewangle); xlim([-plim/2 plim/2]); ylim([-plim/2 plim/2]); zlim([-plim/2 plim/2]);
        %drawnow
        %pause(0.05);

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

        h = figure(1);
        set(h,'Visible','on'); 
        set(gca,'nextplot','replacechildren');
        %set(gcf,'Renderer','OpenGL');
        zoom=5;
        plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane); 
        frame = getframe;
        writeVideo(vidObj,frame); 

        %mov=addframe(mov,getframe(fig));
        %drawnow

        data=zeros(totalsteps,1);

        totaltime=0;
        dt=min(dt,dt0);
        mdold=10;
        counter=0;

        %%%%% Parameters for Stress Strain Curve

        glideplane=[ 1 0 0 ]/norm([ 1 0 0 ]);
        num_tot_points=10;
        appliedstress_basis = appliedstress; %13.0676e-14.*(1/(2*sqrt(6))).*([2 0 1; 0 -2 -1; 1 -1 0]);
        %stress_step=5;
        starting_stress=0;
        ep_tot=zeros(3,3);
        ep_tot_conv=zeros(3,3);
        norm_shear_strain_old=0;
        norm_shear_strain=0;
        norm_shear_strain_conv=0;
        curstep=1;
        %total_shear=1e-4;
        j=0;
        %for i=1:num_tot_points
        i=0;
        %load restart.mat
        %p_file=fopen(filename,'w');
        %vidObj = VideoWriter(videoname);
        %vidObj.Quality = 100;
        %vidObj.FrameRate = 10;
        %open(vidObj);
        %h = figure(1);
        %set(h,'Visible','on'); 
        %set(gca,'nextplot','replacechildren');
        %set(gcf,'Renderer','OpenGL');
        %plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane); view(viewangle); xlim([-plim/2 plim/2]); ylim([-plim/2 plim/2]); zlim([-plim/2 plim/2]);
        %frame = getframe;
        %writeVideo(vidObj,frame); 
        rn_conv=rn;
        links_conv=links;
        linksnew=links;
        while(norm_shear_strain<total_shear)

            if i==0
                stress_multi=0;
            else
                stress_multi=i*stress_step + starting_stress;
            end
            appliedstress=stress_multi*appliedstress_basis;
            curstep=curstep + 1;
            totaltime=0;
            i=i+1;
            convstep=0;

        %         if i>10
        %             doremesh=0;
        %             ducollision=0;
        %             doseparation=0;
        %         end  
            mini = 1;
            while(1)
                j=j+1;
                convstep=convstep+1;
                %integrating equation of motion
                [rnnew,vn,dt,fn,fseg,totaltime]=feval(integrator,rn,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility,dopartials,stackforcevec,totaltime,rann,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,curstep,printfreq,doshielding,SurfacePlane,dobox,boxD);
                % plastic strain and plastic spin calculations
                [ep_inc,wp_inc,L]=calcplasticstrainincrement(rnnew,rn,linksnew,links,(2*plim)^3);
                [ep_inc_conv,wp_inc_conv,L_conv]=calcplasticstrainincrement(rnnew,rn_conv,linksnew,links_conv,(2*plim)^3);

                plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane); %view(viewangle); xlim([-plim/zoom plim/zoom]); ylim([-plim/zoom plim/zoom]); zlim([-plim/zoom plim/zoom]);
                drawnow 

                [ep_tot,norm_shear_strain,conv,ep_tot_conv,norm_shear_strain_conv,mini]=TotalStrain(ep_inc,glideplane,ep_tot,stress_multi,norm_shear_strain_old, totaltime,total_shear,plim,j,p_file,ep_inc_conv,wp_inc_conv,ep_tot_conv,mini,dt,eps,L_conv,convstep);
                norm_shear_strain_old=norm_shear_strain;

                if(conv==1)

                    rn_conv=rn;
                    links_conv=links;

                    disp(sprintf('step%3d dt=%e\ttime=%e',...
                        curstep,dt,totaltime));
                    for h=1:numprint
                        if(printnode(h)>vnmax)
                            printnode(h)=vnmax-(h-1);
                        end
                        disp(sprintf('v%d=(%e,%e,%e)\nforce%d=(%e, %e, %e)\n',...
                            printnode(h),vn(printnode(h),1),vn(printnode(h),2),vn(printnode(h),3),printnode(h),fn(printnode(h),1),fn(printnode(h),2),fn(printnode(h),3)));
                    end
                    %fig = figure(); 
                    %plim=max(max(abs(rn(:,1:3))))+20;
                    plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane); %view(viewangle); xlim([-plim/zoom plim/zoom]); ylim([-plim/zoom plim/zoom]); zlim([-plim/zoom plim/zoom]);
                    frame = getframe;
                    writeVideo(vidObj,frame); 
                    %mov=addframe(mov,getframe(fig));
                    %plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane); view(viewangle); xlim([-plim/5 plim/5]); ylim([-plim/5 plim/5]); zlim([-plim/5 plim/5]);
                    %drawnow
                    %pause(0.05);
                    %aviobj=addframe(aviobj,getframe);
                end

                if(conv==1)
                    break;
                end
                numprint=length(printnode);
                vnmax=size(vn,1);

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
                    [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=separationnew(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,mobility,MU,NU,a,Ec,rann,appliedstress,dopartials,stackforcevecnew,rann,gamma,rntol,maxconnections,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD);

                end
                if(docollision)
                    %collision detection and handling
                    [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=collision(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,appliedstress,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD);
                end

                if(doremesh)
                    %remesh
                    [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=remesh(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,appliedstress,mobility,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD);      
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
                %save restart
                consistencycheck(rn,links,connectivity,linksinconnect,rntol);
                conservstakforce(stackforcevec,links,rn,connectivity);

            end
        end
        plotnodes(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane); view(viewangle); xlim([-plim/5 plim/5]); ylim([-plim/5 plim/5]); zlim([-plim/5 plim/5]);
        frame = getframe;
        writeVideo(vidObj,frame); 
        fclose(p_file);
        close(vidObj);
        save restart
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%