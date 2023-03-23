%This function analise the possiblity of the reaction partial=partial+stair-rod and performs it if needed

function [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew,crossslip]=reactionpartialstairod(rn,links,connectivity,linksinconnect,fseg,stackforcevec,burgv,tangent,nodeid,linkid,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,crossslip,integrator,dt,dt0,rmax,totaltime,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding)


eps=1e-6;
maxconnections=8;%Can take it from the input
delta=1e-16;
numtang=size(tangent,1);
%rntol=1e-2;
%atol=1e-3;

normstairrodacute=norm(gamma.*(1/sqrt(3))*[2 2 0]);
Ethresholdacute=2.13e1;%0.4e1;%1.8e1; %J/a
normstairrodobtuse=norm(gamma.*(1/sqrt(3))*[2 0 0]);
Ethresholdobtuse=1e20; %J/a
PowerNucleation=0;
rnold=rn;
linksold=links;
connectivityold=connectivity;
linksinconnectold=linksinconnect;
fsegold=fseg;
stackforcevecold=stackforcevec;
rnnew=rn;
linksnew=links;
connectivitynew=connectivity;
linksinconnectnew=linksinconnect;
fsegnew=fseg;
stackforcevecnew=stackforcevec;
c=connectivity(nodeid,1);
perfmatrix=[1 -1 0;-1 0 1;0 -1 1;0 1 1;1 0 1;1 1 0]./norm([1 1 0]);%unit perfect vector matrix
 % calculate the power dissipated by the connected geometry and not splitting
        
for i=1:numtang
    glideplane(i,:)=cross(burgv(i,:),tangent(i,:));
    if(norm(glideplane(i,:))>0)
        glideplane(i,:)=glideplane(i,:)/norm(glideplane(i,:));
    else
        glideplane(i,:)=stackforcevec(linkid(i),:)/norm(stackforcevec(linkid(i),:));
    end
    posnodeinsegment(i)=connectivity(nodeid,2*i+1);%This is the position of the node in the segment, can be 1 or 2
    othernodeid(i)=links(linkid(i),3-posnodeinsegment(i));
end
nodelist=[nodeid othernodeid]';
numnodes=length(nodelist);
cmax=size(connectivity,2);
clist=zeros(numnodes,cmax);
clist(1:numnodes,1)=[connectivity(nodelist,1)];
for k=1:numnodes
    clist(k,2:1+connectivity(nodelist(k),1))= linspace(1,connectivity(nodelist(k),1),connectivity(nodelist(k),1));
end
[vntmp,fntmp]=feval(mobility,fseg,rn,links,connectivity,nodelist,clist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad); 
Powermax=0;
for k=1:numnodes
    Powermax=Powermax+dot(vntmp(k,:),fntmp(k,:));
end
numlinks=length(linkid);

if(norm(cross(burgv(1,:),burgv(2,:)))<eps)%
    burgv=burgv(1,:);
    if((norm(glideplane(1,:))>eps)&(norm(glideplane(2,:))>eps)&(norm(glideplane(1,:)-glideplane(2,:))<eps))%If not for sure there is not going to be cross-slip
        
        numperfvec=size(perfmatrix,1);
        for i=1:numperfvec
            for j=1:numtang
                distoperf(j)=norm(cross(tangent(j,:)/norm(tangent(j,:)),perfmatrix(i,:)));
            end
            match=0;
            for j=1:numtang
                if(distoperf(j)<=rntol)
                    match=match+1;
                end
            end
            if((match==numtang)&(cross(tangent(1,:)/norm(tangent(1,:)),tangent(2,:)/norm(tangent(2,:)))<=rntol))%So both arms are close enough and in the same direction
                for j=1:numtang
                    if(posnodeinsegment(j)==2)
                        tangent(j,:)=-tangent(j,:);
                    end
                    rn(othernodeid(j),1:3)=(dot(tangent(j,:),perfmatrix(i,:))).*perfmatrix(i,:)+rn(nodeid,1:3);%Now I have the nodes in a perfect line
                    newtangunit(j,:)=(rn(links(linkid(j),2),1:3)-rn(links(linkid(j),1),1:3))/norm(rn(links(linkid(j),2),1:3)-rn(links(linkid(j),1),1:3));%Here I have the real direction of the line
                end
                rnnewpos=rn;
                [newburpar,newstairod,splitdir]=tablepartialstairod(burgv,newtangunit(1,:),rntol);
                if((norm(newburpar)~=0)&(norm(newstairod)~=0)&(norm(splitdir)~=0))
                    newglideplane=cross(splitdir,newtangunit(1,:))/norm(cross(splitdir,newtangunit(1,:)));
                    stairodplane=cross(newstairod,newtangunit(1,:))/norm(cross(newstairod,newtangunit(1,:)));
                    [lrn,lrn2]=size(rn);
                    newnodeid=lrn+1;
                    rnnewpos(newnodeid,:)=zeros(1,lrn2);
                    rnnewpos(newnodeid,1:3)=rn(nodeid,1:3);
                    linklist=[];
                    for j=1:numtang
                        [numlinks,numcolinks]=size(links);
                        newlinkid(j)=numlinks+1;
                        links(newlinkid(j),:)=zeros(1,numcolinks);
                        stackforcevec(newlinkid(j),:)=zeros(1,3);
                        if(posnodeinsegment(j)==2)
                            links(newlinkid(j),1)=othernodeid(j);
                            links(newlinkid(j),2)=newnodeid;
                        else
                            links(newlinkid(j),1)=newnodeid;
                            links(newlinkid(j),2)=othernodeid(j);
                        end
                        links(newlinkid(j),3:5)=newburpar;
                        links(newlinkid(j),6:8)=newglideplane;
                        links(linkid(j),3:5)=newstairod;
                        links(linkid(j),6:8)=stairodplane;
                        linklist=[linklist newlinkid(j) linkid(j)];
                    
                    end
                    linksnewpos=links;
                    splitdirmatrix=[splitdir;-splitdir];
                    numsplitdir=size(splitdirmatrix,1);
                    numnewlink=length(newlinkid);
                    for j=1:numnewlink
                        fseg=[fseg;zeros(1,6)];
                    end
                    stackforcenewpos=stackforcevec;
                    for j=1:numsplitdir
                        PowerNucleation=0;
                        rn=rnnewpos;
                        links=linksnewpos;
                        stackforcevec=stackforcenewpos;
                        rn(newnodeid,1:3)=rn(newnodeid,1:3)+(2*rann).*splitdirmatrix(j,:);%New position of the node
                        newtangpart=zeros(2,3);
                        [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);                    
                    
                        for k=1:c
                            newtangpart(k,:)=(rn(links(newlinkid(k),2),1:3)-rn(links(newlinkid(k),1),1:3))/norm((rn(links(newlinkid(k),2),1:3)-rn(links(newlinkid(k),1),1:3)));
                            stackforcevec(newlinkid(k),:)=gamma.*cross(newtangpart(k,:),splitdirmatrix(j,:))/norm(cross(newtangpart(k,:),splitdirmatrix(j,:)));
                            gammastairod(k,:)=stackforcevec(linkid(k),:)-stackforcevec(newlinkid(k),:);
                            stackforcevec(linkid(k),:)=gammastairod(k,:);
                            if(abs(norm(stackforcevec(linkid(k),:))-normstairrodacute)<=delta)
                                L=norm(rn(links(linkid(k),2),1:3)-rn(links(linkid(k),1),1:3));
                                PowerNucleation=PowerNucleation+Ethresholdacute*L;
                            elseif(abs(norm(stackforcevec(linkid(k),:))-normstairrodobtuse)<=delta)
                                L=norm(rn(links(linkid(k),2),1:3)-rn(links(linkid(k),1),1:3));
                                PowerNucleation=PowerNucleation+Ethresholdobtuse*L;
                            end
                        end
                        fseg=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,appliedstress,linklist,mobility,dopartials,stackforcevec,doinclusion,inclusion_pos_rad);             
                        nodelist=[nodeid newnodeid othernodeid]';
                        numnodes=length(nodelist);
                        cmax=size(connectivity,2);
                        clist=zeros(numnodes,cmax);
                        clist(1:numnodes,1)=[connectivity(nodelist,1)];
                        for k=1:numnodes
                            clist(k,2:1+connectivity(nodelist(k),1))= linspace(1,connectivity(nodelist(k),1),connectivity(nodelist(k),1));
                        end
                        [vntmp,fntmp]=feval(mobility,fseg,rn,links,connectivity,nodelist,clist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad); 
                        for k=1:numnodes
                            rn(nodelist(k),4:6)=vntmp(k,:);
                        end
                                                          
                    %Check the area change
                    
                        vec1=rn(newnodeid,1:3)-rn(othernodeid(1),1:3);
                        vec2=rn(newnodeid,1:3)-rn(othernodeid(2),1:3);
                        proj1=dot(rn(newnodeid,4:6),vec1);
                        proj2=dot(rn(newnodeid,4:6),vec2);
                        vec3=(rn(othernodeid(1),1:3)-rn(othernodeid(2),1:3))/norm(rn(othernodeid(1),1:3)-rn(othernodeid(2),1:3));
                        vec4=(rn(othernodeid(1),1:3)-rn(newnodeid,1:3))/norm(rn(othernodeid(1),1:3)-rn(newnodeid,1:3));
                        proj3=dot(rn(othernodeid(1),4:6),vec3);
                        proj4=dot(rn(othernodeid(1),4:6),vec4);
                        vec5=(rn(othernodeid(2),1:3)-rn(othernodeid(1),1:3))/norm(rn(othernodeid(2),1:3)-rn(othernodeid(1),1:3));
                        vec6=(rn(othernodeid(2),1:3)-rn(newnodeid,1:3))/norm(rn(othernodeid(2),1:3)-rn(newnodeid,1:3));
                        proj5=dot(rn(othernodeid(2),4:6),vec5);
                        proj6=dot(rn(othernodeid(2),4:6),vec6);
                        darea2dt=proj1+proj2+proj3+proj4+proj5+proj6;
                                                          
                    %Spacing increment
                    
                     %   [lrn,lrn2]=size(rn);
                    %dist1=rn(newnodeid,1:3)-rn(nodeid,1:3);
                    %integrator='int_eulerbackward';
                    %[rnnewtmp,vn,dt,fn,fsegnewtmp]=feval(integrator,rn(:,[1 2 3 lrn2]),dt,dt0,MU,NU,a,Ec,links,connectivity,appliedstress,rmax,rntol,mobility,dopartials,stackforcevec,totaltime);
                    %dist2=rnnewtmp(newnodeid,1:3)-rnnewtmp(nodeid,1:3);
                        Powertest=0;
                        for k=1:numnodes
                            Powertest=Powertest+dot(vntmp(k,:),fntmp(k,:));
                        end 
                        if((Powertest>(Powermax+PowerNucleation))&(darea2dt>0))   %&(dot(dist1,dist2)>0)&((norm(dist2)-norm(dist1))>eps))
                        
                            crossslip=1;
                            rnnew=[];
                            linksnew=[];
                            stackforcevecnew=[];
                            fsegnew=[];
                            connectivitynew=[];
                            linksinconnectnew=[];
                            rnnew=rn;
                            linksnew=links;
                            stackforcevecnew=stackforcevec;
                            fsegnew=fseg;
                            [connectivitynew,linksinconnectnew]=genconnectivity(rn,links,maxconnections);
                            Powermax=Powertest-PowerNucleation;
                            disp('Fleisher is taking place');
                        
                        end
                    end
                    rn=rnold;
                    links=linksold;
                    stackforcevec=stackforcevecold;
                    fseg=fsegold;
                    [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
                                           
                end          
            end
        end
    end
end
