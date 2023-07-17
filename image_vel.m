function [ vn ] = image_vel( rn,vn,SurfacePlane )
%Velocity projected on to a surface plane
%   If the flag of one node is 4 that means that the node is at a surface
%   and then itis constraint to move along the surface plane. The nornal

    nnodes=length(rn(:,4));
    SurfacePlane = SurfacePlane(1:3);
    for i=1:nnodes,
        flag = rn(i,4);
        if( flag==-4 || flag==4)
            nsurf=norm(SurfacePlane);
            UnitSurf=SurfacePlane./nsurf;
            vn(i,:)=-cross(cross(vn(i,:),UnitSurf),UnitSurf);
%         elseif flag==4
%             nsurf=norm(SurfacePlane);
%             UnitSurf=-SurfacePlane./nsurf;
%             vn(i,:)=-cross(cross(vn(i,:),UnitSurf),UnitSurf);
        end
    end

end

