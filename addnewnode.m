%This function add a new node for the general splitting function

function [rn,links,stackforcevec,newnodeid,newlinkid]=addnewnode(rn,links,stackforcevec,tangent,linkid,nodeid,posnodeinsegment,rann,dist)

if(dist==[])
    dist=2*rann;
end

[numnodes,flagpos]=size(rn);
newnodeid=numnodes+1;
rn=[rn;zeros(1,flagpos)];
posnewnode=rn(nodeid,1:3)+dist.*tangent/norm(tangent);
rn(newnodeid,1:3)=posnewnode;
[numlinks,columlinks]=size(links);
newlinkid=numlinks+1;
links=[links;zeros(1,columlinks)];
if(posnodeinsegment==1)
    links(newlinkid,1)=nodeid;
    links(newlinkid,2)=newnodeid;
    links(linkid,1)=newnodeid;
elseif(posnodeinsegment==2)
    links(newlinkid,1)=newnodeid;
    links(newlinkid,2)=nodeid;
    links(linkid,2)=newnodeid;
end
links(newlinkid,3:5)=links(linkid,3:5);%Same Burgers vector
links(newlinkid,6:8)=links(linkid,6:8);%Same glide plane
stackforcevec=[stackforcevec;zeros(1,3)];
stackforcevec(newlinkid,:)=stackforcevec(linkid,:);%Same stackvector


                            

