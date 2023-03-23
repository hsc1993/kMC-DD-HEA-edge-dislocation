function [rn,vn,dt,fn,fseg,DDtime]=int_eulerbackward(rn,dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility,dopartials,stackforcevec,DDtime,rann,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,curstep,printfreq,doshielding,SurfacePlane,dobox,boxD,Bt_edge,Bt_screw,edgecrss,screwcrss)

%implicit numerical integrator using the Backward Euler or Trapezoid method
%dt: suggested timestep (usually from previous iteration)
%dt0: maximum allowed timestep

%dummy variable
t=0;
rnold=rn;

%scramble rn into a single column vector
rnvec=[rn(:,1);rn(:,2);rn(:,3)];
flag=rn(:,4);

%Backward Euler
rnvec0=rnvec;


[vnvec0,fn,fseg,flag]=drndt(t,rnvec0,flag,MU,NU,a,Ec,links,connectivity,appliedstress,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD,dt,Bt_edge,Bt_screw,edgecrss,screwcrss);

%dt=1/max(1/dt0,max(vnvec0)/rmax);
%dt=dt0;

dt1=dt;
maxiter=1;
convergent=0;
rnvec1=[];
while(~convergent)
    rnvec1=rnvec0+(vnvec0.*dt);
    for iter=1:maxiter
        
        [vnvec,fn,fseg,flag]=drndt(t,rnvec1,flag,MU,NU,a,Ec,links,connectivity,appliedstress,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD,dt,Bt_edge,Bt_screw,edgecrss,screwcrss);
        
        err=rnvec1-rnvec0-vnvec.*dt;          %backward Euler
        %err=rnvec1-rnvec0-(vnvec+vnvec0)./2.*dt; %trapzoid
%         errmag=max(sqrt(err(:,1).^2+err(:,2).^2+err(:,3).^2));
        errmag=max(abs(err));
        jump=rnvec1-rnvec0;
        %jumpmag=max(sqrt(jump(:,1).^2+jump(:,2).^2+jump(:,3).^2));
        jumpmag=max(abs(jump));
        %if(mod(curstep,printfreq)==0)
        %disp(sprintf('iter=%d err=%e jump=%e',iter,errmag,jumpmag));
        
        %end
        if((errmag<rntol))%&&(jumpmag<=(rann))
            convergent=1;
            break;
        else
            rnvec1=rnvec1-err;
        end
    end
    %disp(sprintf('ppt=%i %i\n',size(ppt,1),size(ppt,2)));
    if(convergent)
        break;
    else
        dt=dt/2;
    end
end

%unscramble rn and vn vectors
rn=[reshape(rnvec1,length(rnvec1)/3,3),flag];
vn=reshape(vnvec,length(vnvec)/3,3); % trapezoidal rule modificatio
vn0=reshape(vnvec0,length(vnvec)/3,3); % trapezoidal rule modificatio

errvec=reshape(err,length(err)/3,3);    

num_inclusions=size(inclusion_pos_rad,1);
for i=1:num_inclusions
   rn(:,4)=rn(:,4)+(rn(:,4)==0 & (sqrt(sum((rn(:,1:3)-repmat(inclusion_pos_rad(i,1:3),size(rn,1),1)).^2,2))<repmat(inclusion_pos_rad(i,4),size(rn,1),1)))*3;
end

% for n=1:size(rn(:,1),1)
%     if norm(rn(n,1:3)-inclusion_pos_rad(1,1:3))<=inclusion_pos_rad(1,4) || rn(n,4)==3
%         %rn(n,4)=3;
%         disp(sprintf('\nnode=%i flag = %i rad=%d dist=%d',n,rn(n,4),inclusion_pos_rad(1,4),norm(rn(n,1:3)-inclusion_pos_rad(1,1:3))));
%         disp(sprintf('vn=%d %d %d \n rn=%d %d %d\n',vn(n,1),vn(n,2),vn(n,3),rn(n,1),rn(n,2),rn(n,3)));
%         pause(0.2)
%     end
% end
%automatically adjust time step
if((dt==dt1)&(iter==1))
    factor=1.2*(1/(1+(1.2^20-1)*(errmag/rntol)))^0.05;
    dt=min(dt1*factor,dt0);
end
%disp(sprintf('totaltime=%e dt=%e\n',totaltime,dt));

DDtime=DDtime+dt;
% totaltime=totaltime+dt;
