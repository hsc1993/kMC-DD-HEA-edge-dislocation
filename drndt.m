function [vnvec,fn,fseg,flag] = drndt(t,rnvec,flag,MU,NU,a,Ec,links,connectivity,appliedstress,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD,dt,Bt_edge,Bt_screw,edgecrss,screwcrss)


%unscramble rn
rn=[reshape(rnvec,length(rnvec)/3,3),flag];

%rn(:,1:3)

%nodal driving force
fseg=segforcevec(MU,NU,a,Ec,rn,links,appliedstress,[],mobility,dopartials,stackforcevec,doinclusion,inclusion_pos_rad,doSFT,SFT_plane,doshielding,SurfacePlane,dobox,boxD);
% fseg [appliedstress (N/A^2)]*[A]*[A] -> N

% Alan edits to introduce crss 
b = [2.7973 0 0];
s = b/norm(b);
theta = zeros(size(rn,1)-1,1);
n = [0 1 0];

for ii = 1:size(rn,1)-1
    l = rn(ii+1,1:3)-rn(ii,1:3);
    t = l/norm(l);
    theta(ii) = real(acos(b*t'));
    crss_norm = abs(edgecrss*sin(theta(ii))+screwcrss*cos(theta(ii)));  %[N/A^2]
    delta_f = crss_norm*norm(b)*norm(l);
    crss_direction = cross(n,t);
    crss_vector = crss_norm*crss_direction;

    % use projection of crss on each dimension to calculate the available
    % force in each dimension (most of the time its just along [100] direction)
    for dim = 1:3
        if fseg(ii,dim)>0
            fseg(ii,dim) = fseg(ii,dim) - crss_vector(dim);
                if fseg(ii,dim)<0
                    fseg(ii,dim) = 0;
                end
        end
        if fseg(ii,dim)<0
            fseg(ii,dim) = fseg(ii,dim) + crss_vector(dim);
                if fseg(ii,dim)>0
                    fseg(ii,dim) = 0;
                end
        end
        if fseg(ii,dim+3)>0
            fseg(ii,dim+3) = fseg(ii,dim+3) - crss_vector(dim);
                if fseg(ii,dim+3)<0
                    fseg(ii,dim+3) = 0;
                end
        end
        if fseg(ii,dim+3)<0
            fseg(ii,dim+3) = fseg(ii,dim+3) + crss_vector(dim);
                if fseg(ii,dim+3)>0
                    fseg(ii,dim+3) = 0;
                end
        end
    end
end

% crss

%mobility function
% [vn,fn]=feval(mobility,fseg,rn,links,connectivity,[],[],mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,dt);
[vn,fn]=feval(mobility,fseg,rn,links,connectivity,[],[],mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,dt,Bt_edge,Bt_screw); % Alan uses this to introduce Bt dependence on Temperature
% [vn,fn]=feval(mobility,fseg,rn,links,connectivity,[],[]);  % Alan use this to study bcc edge model

%fixed nodes (flag~=0) are not allowed to move
vn=vn.*((rn(:,4)==0 | rn(:,4)==3)*[1 1 1]);
%nodes on precipitate are not allowed to move
% num_inclusions=size(inclusion_pos_rad,1);
% for i=1:num_inclusions
%     vn=vn.*((sqrt(sum((rn(:,1:3)-repmat(inclusion_pos_rad(i,1:3),size(rn,1),1)).^2,2))>repmat(inclusion_pos_rad(i,4),size(rn,1),1))*[1 1 1]);
% end
    
%if flag==4 then is a surface node and can move only in the surface plane.
%The normal to the surface plane is estimated with the links to the surface
%node
vn = image_vel(rn,vn,SurfacePlane);
flag=rn(:,4);
%make up the last
vnvec=reshape(vn,length(vn(:,1))*3,1);

%vn
%disp(sprintf('t=%20.10e vn(1,:)=(%20.10e,%20.10e,%20.10e)',t,vn(1,1),vn(1,2),vn(1,3)));
%pause(0.2)
%if(t>0)
%   pause
%end