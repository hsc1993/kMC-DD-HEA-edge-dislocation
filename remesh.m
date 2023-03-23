function [rn,links,connectivity,linksinconnect,fseg,stackforcevec]=remesh(rn,links,connectivity,linksinconnect,fseg,lmin,lmax,areamin,areamax,MU,NU,a,Ec,appliedstress,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD,dt,Bt_edge,Bt_screw)
dt = 0;  % dummy 
% first coarsen the parts of the mesh than can be coarsened
% then refine the parts of the mesh that need it
% meshcoarsen is done first because it relies on positions of nodes that were introduced at the end of the previous time step
% do not change the order of these two subroutines
[rn,links,connectivity,linksinconnect,fseg,stackforcevec]=meshcoarsen(rn,links,connectivity,linksinconnect,fseg,lmin,lmax,areamin,MU,NU,a,Ec,appliedstress,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD,dt,Bt_edge,Bt_screw);
[rn,links,connectivity,linksinconnect,fseg,stackforcevec]=meshrefine(rn,links,connectivity,linksinconnect,fseg,lmin,lmax,areamax,MU,NU,a,Ec,appliedstress,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD,dt,Bt_edge,Bt_screw);
[rn]=surface(rn,SurfacePlane);


function [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=meshcoarsen(rn,links,connectivity,linksinconnect,fseg,lmin,lmax,areamin,MU,NU,a,Ec,appliedstress,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD,dt,Bt_edge,Bt_screw)
rnnew=rn;
[lrn lrn2]=size(rn);
linksnew=links;
connectivitynew=connectivity;
linksinconnectnew=linksinconnect;
fsegnew=fseg;
stackforcevecnew=stackforcevec;
normstairrod=norm(1/6.*[1 1 0]);
normperfect=norm(1/2.*[1 1 0]);
areamin2=areamin*areamin;
i=1;
eps=1e-12;
delta=1e-16;
lmin_orig=lmin;
lmax_orig=lmax;
while i<=length(rnnew(:,1))
    if (connectivitynew(i,1)==2 && rnnew(i,lrn2)==0)
        % the discretization node is normal so set up the conditions to check for whether is should be coarsened away
        link1=connectivitynew(i,2);
        link2=connectivitynew(i,4);
        if((abs(norm(linksnew(link1,3:5))-normstairrod)>eps)&(abs(norm(linksnew(link2,3:5))-normstairrod)>eps))
            posi_link1=connectivitynew(i,3);
            posi_link2=connectivitynew(i,5);
            posnoti_link1=3-posi_link1;
            posnoti_link2=3-posi_link2;
            link1_nodenoti=linksnew(link1,posnoti_link1);
            link2_nodenoti=linksnew(link2,posnoti_link2);
            vec1=rnnew(link1_nodenoti,1:3)-rnnew(i,1:3);
            vec2=rnnew(link2_nodenoti,1:3)-rnnew(i,1:3);
            vec3=vec2-vec1;
            r1=sqrt(vec1*vec1');
            r2=sqrt(vec2*vec2');
            r3=sqrt(vec3*vec3');
%             if doinclusion
%                 if(norm(rnnew(i,1:3)-inclusion_pos_rad(1,1:3))<(inclusion_pos_rad(1,4))+3 || ...
%                         norm(rnnew(link1_nodenoti,1:3)-inclusion_pos_rad(1,1:3))<(inclusion_pos_rad(1,4))+3 || ...
%                         norm(rnnew(link2_nodenoti,1:3)-inclusion_pos_rad(1,1:3))<(inclusion_pos_rad(1,4))+3)
%                    lmin = 2;
%                    lmax = 4;
%                    rann=1.75*(1/sqrt(2));
%                    rntol = rann/2;
%                    areaminmag(1)=2*lmin*rann*cos(asin(2*rann/(2*pi*lmin))); 
%                    areaminmag(2)=2*lmax*rntol*cos(asin(2*rntol/(2*pi*lmax)));
%                    areamin=min(areaminmag);
%                    areamin2=areamin*areamin;
%                 else
%                     lmin=lmin_orig;
%                     lmax=lmax_orig;
%                 end
%             end
            if r3<lmax %if the coarsening would result in a link length larger than lmax don't coarsen
                s=0.5*(r1+r2+r3);
                area2=(s*(s-r1)*(s-r2)*(s-r3));
                dvec1dt=rnnew(link1_nodenoti,4:6)-rnnew(i,4:6);
                dvec2dt=rnnew(link2_nodenoti,4:6)-rnnew(i,4:6);
                dvec3dt=dvec2dt-dvec1dt;
                dr1dt=(vec1*dvec1dt')/(r1+delta);
                dr2dt=(vec2*dvec2dt')/(r2+delta);
                dr3dt=(vec3*dvec3dt')/(r3+delta);
                dsdt=0.5*(dr1dt+dr2dt+dr3dt);
                darea2dt=dsdt*(s-r1)*(s-r2)*(s-r3);
                darea2dt=darea2dt+s*(dsdt-dr1dt)*(s-r2)*(s-r3);
                darea2dt=darea2dt+s*(s-r1)*(dsdt-dr2dt)*(s-r3);
                darea2dt=darea2dt+s*(s-r1)*(s-r2)*(dsdt-dr3dt);
                if ((area2<areamin2)&(darea2dt<0.0d0))|((r1<lmin)|(r2<lmin))%remove single node critierion
                    %the area is less than minimum and shrinking or one of the arms is less than lmin and shrinking
                    [rnnew,connectivitynew,linksnew,linksinconnectnew,fsegnew,mergednodeid,stackforcevecnew]=mergenodes(rnnew,connectivitynew,linksnew,linksinconnectnew,fsegnew,link2_nodenoti,i,MU,NU,a,Ec,mobility,dopartials,stackforcevecnew,inclusion_pos_rad);
                    if mergednodeid>0
                        for j=1:connectivitynew(mergednodeid,1)
                            linkm=connectivitynew(mergednodeid,2*j);
                            posm=connectivitynew(mergednodeid,2*j+1);
                            connode=linksnew(linkm,3-posm);
                            if ((connode==i)|(connode==link1_nodenoti))
                                fsegnew=segforcevec(MU,NU,a,Ec,rnnew(:,[1 2 3 lrn2]),linksnew,appliedstress,linkm,mobility,dopartials,stackforcevecnew,doinclusion,inclusion_pos_rad,doSFT,SFT_plane,doshielding,SurfacePlane,dobox,boxD);
                                for j=1:2
                                    nodelist=linksnew(linkm,j);
                                    clist=[connectivitynew(nodelist,1) linspace(1,connectivitynew(nodelist,1),connectivitynew(nodelist,1))];
                                    [rnnew(nodelist,4:6),fntmp]=feval(mobility,fsegnew,rnnew,linksnew,connectivitynew,nodelist,clist,mobility,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,dt,Bt_edge,Bt_screw);
                                end    
                            end
                        end
                    end
                else %((area2>=areamin2)|(dareadt>=0.0d0)) & ((r1>=lmin)|(dr1dt>=lmin)) & ((r2>=lmin)|(dr2dt>=lmin))
                i=i+1;
                end
            else % r3>=lmax
                i=i+1;
            end
        else
        i=i+1;    
    end
    else % connectivitynew(i,1)>2
        i=i+1;
    end
       
end %while loop




function [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=meshrefine(rn,links,connectivity,linksinconnect,fseg,lmin,lmax,areamax,MU,NU,a,Ec,appliedstress,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD,dt,Bt_edge,Bt_screw)
[lrn lrn2]=size(rn);
lrn3=lrn2-1;
rnnew=rn;
linksnew=links;
connectivitynew=connectivity;
linksinconnectnew=linksinconnect;
fsegnew=fseg;
stackforcevecnew=stackforcevec;
normstairrod=norm(1/6.*[1 1 0]);
normperfect=norm(1/2.*[1 1 0]);
lmin2=2*lmin;
areamax2=areamax*areamax;
lmin_orig=lmin;
lmax_orig=lmax;
for i=1:lrn
%     if doinclusion
%         if(norm(rnnew(i,1:3)-inclusion_pos_rad(1,1:3))<(inclusion_pos_rad(1,4))+6)
%            lmin = 2;
%            lmax = 3;
%         else
%             lmin=lmin_orig;
%             lmax=lmax_orig;
%         end
%     end
    if (connectivitynew(i,1)==2) & ((rnnew(i,lrn2)==0) || (rnnew(i,lrn2)==4) || (rnnew(i,lrn2)==-4)) 
        firstconnection=1;
        secondconnection=2;
        link1=connectivitynew(i,2*firstconnection);
        link2=connectivitynew(i,2*secondconnection);
        
        if((abs(norm(linksnew(link1,3:5))-normstairrod)>eps)&(abs(norm(linksnew(link2,3:5))-normstairrod)>eps))
            posi_link1=connectivitynew(i,3);
            posi_link2=connectivitynew(i,5);
            posnoti_link1=3-posi_link1;
            posnoti_link2=3-posi_link2;
            link1_nodenoti=linksnew(link1,posnoti_link1);
            link2_nodenoti=linksnew(link2,posnoti_link2);
            vec1=rnnew(link1_nodenoti,1:3)-rnnew(i,1:3);
            vec2=rnnew(link2_nodenoti,1:3)-rnnew(i,1:3);
            vec3=vec2-vec1;
            r1=sqrt(vec1*vec1');
            r2=sqrt(vec2*vec2');
            r3=sqrt(vec3*vec3');
            s=0.5*(r1+r2+r3);
            area2=(s*(s-r1)*(s-r2)*(s-r3));
%             if doinclusion
%                 if(norm(rnnew(i,1:3)-inclusion_pos_rad(1,1:3))<(inclusion_pos_rad(1,4))+3 || ...
%                         norm(rnnew(link2_nodenoti,1:3)-inclusion_pos_rad(1,1:3))<(inclusion_pos_rad(1,4))+3)
%                    lmin = 2;
%                    lmax = 4;
%                    lmin2=2*lmin;
%                    areamax=(1/4)*lmax*lmax;
%                    areamax2=areamax*areamax;
%                 else
%                     lmin=lmin_orig;
%                     lmax=lmax_orig;
%                 end
%                 
%             end
            if (((area2>areamax2)&(r2>=lmin2)&(link2_nodenoti<=lrn))|(r2>lmax))
                %conditions necessary to bisect the second link are met 
                posvel=[(rnnew(i,1:lrn3)+rnnew(link2_nodenoti,1:lrn3))./2];
                flag = 0;
                if(rnnew(i,lrn2)==4 || rnnew(i,lrn2)==-4 || ...
                        rnnew(link2_nodenoti,lrn2)==4 || rnnew(link2_nodenoti,lrn2)==-4)
                    flag=4;
                end
                [rnnew,linksnew,connectivitynew,linksinconnectnew,stackforcevecnew]=splitnode(rnnew,linksnew,connectivitynew,linksinconnectnew,i,secondconnection,posvel,stackforcevecnew,flag);
                newnode=length(rnnew(:,1));
                for j=1:connectivitynew(newnode,1)
                    linkid=connectivitynew(newnode,2*j);
                    oldnode=linksnew(linkid,3-connectivitynew(newnode,2*j+1));
                    fsegnew=segforcevec(MU,NU,a,Ec,rnnew(:,[1 2 3 lrn2]),linksnew,appliedstress,linkid,mobility,dopartials,stackforcevecnew,doinclusion,inclusion_pos_rad,doSFT,SFT_plane,doshielding,SurfacePlane,dobox,boxD);
                    clist=[connectivitynew(oldnode,1) linspace(1,connectivitynew(oldnode,1),connectivitynew(oldnode,1))];
                    [rnnew(oldnode,4:6),fntmp]=feval(mobility,fsegnew,rnnew,linksnew,connectivitynew,oldnode,clist,mobility,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,dt,Bt_edge,Bt_screw);
                end
                clist=[connectivitynew(newnode,1) linspace(1,connectivitynew(newnode,1),connectivitynew(newnode,1))];
                [rnnew(newnode,4:6),fntmp]=feval(mobility,fsegnew,rnnew,linksnew,connectivitynew,newnode,clist,mobility,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,dt,Bt_edge,Bt_screw);
            end
%             if doinclusion
%                 if(norm(rnnew(i,1:3)-inclusion_pos_rad(1,1:3))<(inclusion_pos_rad(1,4))+3 || ...
%                         norm(rnnew(link1_nodenoti,1:3)-inclusion_pos_rad(1,1:3))<(inclusion_pos_rad(1,4))+3)
%                    lmin = 2;
%                    lmax = 4;
%                    areamax=(1/4)*lmax*lmax;
%                    areamax2=areamax*areamax;
%                 else
%                     lmin=lmin_orig;
%                     lmax=lmax_orig;
%                 end
%             end
            if (((area2>areamax2)&(r1>=lmin2)&(link1_nodenoti<=lrn))|(r1>lmax))
                %conditions necessary to bisect the first link are met
                posvel=[(rnnew(i,1:lrn3)+rnnew(link1_nodenoti,1:lrn3))./2];
                flag = 0;
                if(rnnew(i,lrn2)==4 || rnnew(i,lrn2)==-4 || ...
                        rnnew(link2_nodenoti,lrn2)==4 || rnnew(link2_nodenoti,lrn2)==-4)
                    flag=4;
                end
                [rnnew,linksnew,connectivitynew,linksinconnectnew,stackforcevecnew]=splitnode(rnnew,linksnew,connectivitynew,linksinconnectnew,i,firstconnection,posvel,stackforcevecnew,flag);
                newnode=length(rnnew(:,1));
                for j=1:connectivitynew(newnode,1)
                    linkid=connectivitynew(newnode,2*j);
                    oldnode=linksnew(linkid,3-connectivitynew(newnode,2*j+1));
                    fsegnew=segforcevec(MU,NU,a,Ec,rnnew(:,[1 2 3 lrn2]),linksnew,appliedstress,linkid,mobility,dopartials,stackforcevecnew,doinclusion,inclusion_pos_rad,doSFT,SFT_plane,doshielding,SurfacePlane,dobox,boxD);
                    clist=[connectivitynew(oldnode,1) linspace(1,connectivitynew(oldnode,1),connectivitynew(oldnode,1))];
                    [rnnew(oldnode,4:6),fntmp]=feval(mobility,fsegnew,rnnew,linksnew,connectivitynew,oldnode,clist,mobility,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,dt,Bt_edge,Bt_screw);
                end
                clist=[connectivitynew(newnode,1) linspace(1,connectivitynew(newnode,1),connectivitynew(newnode,1))];
                [rnnew(newnode,4:6),fntmp]=feval(mobility,fsegnew,rnnew,linksnew,connectivitynew,newnode,clist,mobility,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,dt,Bt_edge,Bt_screw);
            end
        end
    elseif (connectivitynew(i,1)>2) & ((rnnew(i,lrn2)==0) || (rnnew(i,lrn2)==4) || (rnnew(i,lrn2)==-4)) 
        % check to make sure that no link is larger than lmax
        for j=1:connectivitynew(i,1)
            linkid=connectivitynew(i,2*j);
            if(abs(norm(linksnew(linkid,3:5))-normstairrod)>eps)
                posi=connectivitynew(i,2*j+1);
                nodenoti=linksnew(linkid,3-posi);
                vec1=rnnew(nodenoti,1:3)-rnnew(i,1:3);
                r1=sqrt(vec1*vec1');
                if (r1>lmax)
                    posvel=[(rnnew(i,1:lrn3)+rnnew(nodenoti,1:lrn3))./2];
                    flag = 0;
                    if(rnnew(i,lrn2)==4 || rnnew(i,lrn2)==-4 || ...
                            rnnew(nodenoti,lrn2)==4 || rnnew(nodenoti,lrn2)==-4)
                        flag=4;
                    end
                    [rnnew,linksnew,connectivitynew,linksinconnectnew,stackforcevecnew]=splitnode(rnnew,linksnew,connectivitynew,linksinconnectnew,i,j,posvel,stackforcevecnew,flag);
                    newnode=length(rnnew(:,1));
                    for j=1:connectivitynew(newnode,1)
                        linkid=connectivitynew(newnode,2*j);
                        oldnode=linksnew(linkid,3-connectivitynew(newnode,2*j+1));
                        fsegnew=segforcevec(MU,NU,a,Ec,rnnew(:,[1 2 3 lrn2]),linksnew,appliedstress,linkid,mobility,dopartials,stackforcevecnew,doinclusion,inclusion_pos_rad,doSFT,SFT_plane,doshielding,SurfacePlane,dobox,boxD);
                        clist=[connectivitynew(oldnode,1) linspace(1,connectivitynew(oldnode,1),connectivitynew(oldnode,1))];
                        [rnnew(oldnode,4:6),fntmp]=feval(mobility,fsegnew,rnnew,linksnew,connectivitynew,oldnode,clist,mobility,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,dt,Bt_edge,Bt_screw);
                    end
                    clist=[connectivitynew(newnode,1) linspace(1,connectivitynew(newnode,1),connectivitynew(newnode,1))];
                    [rnnew(newnode,4:6),fntmp]=feval(mobility,fsegnew,rnnew,linksnew,connectivitynew,newnode,clist,mobility,dopartials,stackforcevecnew,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,dt,Bt_edge,Bt_screw);
                end
            end
        end
    end
end



% Check if nodes reaches a surface
function [rnnew]=surface(rn,SurfacePlane)
[lrn lrn2]=size(rn);
rnnew=rn;
for i=1:lrn
   if dot(SurfacePlane,[rnnew(i,1:3) 1])< 0
       rnnew(i,lrn2) = 4;
   end
end