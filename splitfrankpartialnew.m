function [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew,printnode]=splitfrankpartialnew(rn,links,connectivity,linksinconnect,fseg,stackforcevec,burgv,tangent,linkid,glideplane,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,areamin,lmin,printnode,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding)

maxconnections=8;
delta=1e-16;
eps=1e-12;
perfmatrix=[1 -1 0;-1 0 1;0 -1 1;0 1 1;1 0 1;1 1 0]./norm([1 1 0]);
numperfdir=size(perfmatrix,1);
crossplane=[1 1 1; -1 1 -1; -1 -1 1; 1 -1 -1]/norm([1 1 1]);
stackforcefrank=stackforcevec(linkid,:);
normstairrodacute=norm(gamma.*(1/sqrt(3))*[2 2 0]);
normstairrodobtuse=norm(gamma.*(1/sqrt(3))*[2 0 0]);
normperf=norm(1/2.*[1 1 0]);
Ethresholdperf=1e7; %Barrier for the perfect-Shockley dissociation
rnold=rn;
linksold=links;
fsegold=fseg;
stackforcevecold=stackforcevec;
rnnew=rn;
linksnew=links;
connectivitynew=connectivity;
linksinconnectnew=linksinconnect;
fsegnew=fseg;
stackforcevecnew=stackforcevec;

node1=rn(links(linkid,1),1:3);
node2=rn(links(linkid,2),1:3);
L=norm(node1-node2);
nodeid1=links(linkid,1);
nodeid2=links(linkid,2);
nodelist=[nodeid1 nodeid2]';
numnodes=length(nodelist);
cmax=size(connectivity,2);
clist=zeros(numnodes,cmax);
clist(1:numnodes,1)=[connectivity(nodelist,1)];
for k=1:numnodes
    %clist(k,1)=connectivity(nodelist(k),1);
    clist(k,2:1+connectivity(nodelist(k),1))= linspace(1,connectivity(nodelist(k),1),connectivity(nodelist(k),1));
end
[vntmp,fntmp]=feval(mobility,fseg,rn,links,connectivity,nodelist,clist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad); 
Powermax=0;
%vnglide=zeros(numnodes,3);
for k=1:numnodes
    Powermax=Powermax+dot(vntmp(k,:),fntmp(k,:));
end
Powermax=1.01*Powermax;%NO IDEA ABOUT THE BARRIER

[burg1,burg2,splitdir,crossplane]=tablefrankpartialsnew(burgv,tangent/norm(tangent),rntol,stackforcevec(linkid,:));
        
[lrn,lrn2]=size(rn);
newnodeid=lrn+1;
rn(newnodeid,:)=zeros(1,lrn2);
rn(newnodeid,1:3)=(node1+node2)/2;
for j=1:2
    [numlinks,numcolinks]=size(links);
    newlinkid(j)=numlinks+1;
    links(newlinkid(j),:)=zeros(1,numcolinks);
    stackforcevec(newlinkid(j),:)=zeros(1,3);
    if(j==1)
        links(newlinkid(j),1)=nodeid1;
        links(newlinkid(j),2)=newnodeid;
    else
        links(newlinkid(j),1)=newnodeid;
        links(newlinkid(j),2)=nodeid2;
    end
            %links(newlinkid(j),3:5)=burgpart;
            %links(newlinkid(j),6:8)=crossplane;
end
        %links(linkid,3:5)=burgstairrod;
        %links(linkid,6:8)=crossplane;
rnnewpos=rn;
linksnewpos=links;
        %splitdirmatrix=[splitdir;-splitdir];
        %burg1matrix=[burg1;burg1];
        %burg2matrix=[burg2;burg2];
        %crossplanematrix=[crossplane;crossplane];
numsplitdir=size(splitdir,1);
numnewlink=length(newlinkid);
for j=1:numnewlink
    fseg=[fseg;zeros(1,6)];
end
stackforcenewpos=stackforcevec;

check_dir=0;
for i=1:numperfdir
    if(norm(cross(tangent/norm(tangent),perfmatrix(i,:)))<rntol)&(abs(dot(burgv/norm(burgv),tangent/norm(tangent)))<rntol)
        check_dir=1;
    end
end
        
for j=1:numsplitdir
    if((norm(burg1(j,:))>eps)&(norm(burg2(j,:))>eps)&(norm(splitdir(j,:))>eps)&(norm(crossplane(j,:))>eps))
        Ethreshold=0;
        rn=rnnewpos;
        links=linksnewpos;
        stackforcevec=stackforcenewpos;
        rn(newnodeid,1:3)=rn(newnodeid,1:3)+(3*rann).*splitdir(j,:);%New position of the node
        for k=1:2
            links(newlinkid(k),3:5)=burg1(j,:);
            links(newlinkid(k),6:8)=crossplane(j,:);
        end
        links(linkid,3:5)=burg2(j,:);
        links(linkid,6:8)=cross(tangent,burg2(j,:))/norm(cross(tangent,burg2(j,:)));
        newtangpart=zeros(2,3);
        [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);                    
            %stackforcevec(linkid,:)=gamma.*cross(tangent(k,:),(-splitdirmatrix(j,:)))/norm(cross(tangent(k,:),(-splitdirmatrix(j,:))));        
        for k=1:2  %It is not a physical node
            newtangpart(k,:)=(rn(links(newlinkid(k),2),1:3)-rn(links(newlinkid(k),1),1:3))/norm((rn(links(newlinkid(k),2),1:3)-rn(links(newlinkid(k),1),1:3)));
            if(abs(norm(burg2(j,:))-normperf)<eps)
                stackforcevec(newlinkid(k),:)=zeros(1,3);
                Ethreshold=Ethresholdperf;
            else
                stackforcevec(newlinkid(k),:)=gamma.*cross(newtangpart(k,:),splitdir(j,:))/norm(cross(newtangpart(k,:),splitdir(j,:)));
            end 
        end
        stackpart=stackforcevec(newlinkid(k),:);
        stackforcevec(linkid,:)=stackforcevec(linkid,:)-stackpart;
        if(abs(norm(stackforcevec(linkid,:))-normstairrodobtuse)<eps)
            continue;
        end
        linklist=[linkid newlinkid];
        fseg=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,appliedstress,linklist,mobility,dopartials,stackforcevec,doinclusion,inclusion_pos_rad,doSFT,SFT_plane,doshielding);             
       % fseg=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,appliedstress,newlinkid,mobility,dopartials,stackforcevec,doinclusion,inclusion_pos_rad);             
                 
                    
                    %[rn,connectivity,links,linksinconnect,fseg,nodeidmerge,stackforcevec]=mergenodes(rn,connectivity,links,linksinconnect,fseg,othernodeid(1),nodeid,MU,NU,a,Ec,mobility,dopartials,stackforcevec);
        nodelist=[nodeid1 nodeid2 newnodeid]';
        numnodes=length(nodelist);
        cmax=size(connectivity,2);
        clist=zeros(numnodes,cmax);
        clist(1:numnodes,1)=[connectivity(nodelist,1)];
        for k=1:numnodes
            clist(k,1)=connectivity(nodelist(k),1);
            clist(k,2:1+connectivity(nodelist(k),1))= linspace(1,connectivity(nodelist(k),1),connectivity(nodelist(k),1));
        end
                    
        [vntmp,fntmp]=feval(mobility,fseg,rn,links,connectivity,nodelist,clist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad); 
            
        for m=1:numnodes
            rn(nodelist(m),4:6)=vntmp(m,:);
        end
        vec1=(rn(newnodeid,1:3)-rn(nodeid1,1:3))/norm(rn(newnodeid,1:3)-rn(nodeid1,1:3));
        vec2=(rn(newnodeid,1:3)-rn(nodeid2,1:3))/norm(rn(newnodeid,1:3)-rn(nodeid2,1:3));
        proj1=dot(vec1,rn(newnodeid,4:6));
        proj2=dot(vec2,rn(newnodeid,4:6));
        vec3=(rn(nodeid1,1:3)-rn(nodeid2,1:3))/norm(rn(nodeid1,1:3)-rn(nodeid2,1:3));
        vec4=(rn(nodeid1,1:3)-rn(newnodeid,1:3))/norm(rn(nodeid1,1:3)-rn(newnodeid,1:3));
        proj3=dot(rn(nodeid1,4:6),vec3);
        proj4=dot(rn(nodeid1,4:6),vec4);
        vec5=(rn(nodeid2,1:3)-rn(nodeid1,1:3))/norm(rn(nodeid2,1:3)-rn(nodeid1,1:3));
        vec6=(rn(nodeid2,1:3)-rn(newnodeid,1:3))/norm(rn(nodeid2,1:3)-rn(newnodeid,1:3));
        proj5=dot(rn(nodeid2,4:6),vec5);
        proj6=dot(rn(nodeid2,4:6),vec6);
        darea2dt=proj1+proj2;%+proj3+proj4+proj5+proj6;
        Powertest=0;       
                                      
       % if(dot(vntmp(3,1:3),splitdir(j,:))>0)%((proj1>0)|(proj2>0))
            for k=1:numnodes
                Powertest=Powertest+dot(vntmp(k,:),fntmp(k,:));
                    
            end
      %  end
        if(Powertest>(Powermax+Ethreshold))
    
            rnnew=rn;
            linksnew=links;
            fsegnew=fseg;
            stackforcevecnew=stackforcevec;
            [connectivitynew,linksinconnectnew]=genconnectivity(rn,links,maxconnections);
            Powermax=Powertest-Ethreshold;
            rn=rnold;
            links=linksold;
            fseg=fsegold;
            stackforcevec=stackforcevecold;
            disp('Frank partial dissociating');
        end
        rn=rnold;
        links=linksold;
        fseg=fsegold;
        stackforcevec=stackforcevecold; 
                        
    end
end