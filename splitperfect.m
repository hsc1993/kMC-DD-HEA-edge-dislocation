%This function treat the split in partials and cross-slip of a perfect dislocation

function [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew,reset]=splitperfect(rn,links,connectivity,linksinconnect,fseg,stackforcevec,burgv,tangent,nodeid,linkid,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding)

eps=1e-6;
maxconnections=8;%Can take it from the input
delta=1e-16;
rntol=0.5;
numtang=size(tangent,1);
rnnew=rn;
linksnew=links;
connectivitynew=connectivity;
linksinconnectnew=linksinconnect;
fsegnew=fseg;
stackforcevecnew=stackforcevec;
rnold=rn;
linksold=links;
fsegold=fseg;
stackforcevecold=stackforcevec;
[numnodes,flagpos]=size(rn);%Total number of nodes
[numlinks,columlinks]=size(links);%Total number of links
[numfseg,columfseg]=size(fseg);
[numstack,columstack]=size(stackforcevec);

Powermax=0;
reset=0;
c=connectivity(nodeid,1);

if(rann==0)
    rann=a;
end


perfmatrix=[1 -1 0;-1 0 1;0 -1 1;0 1 1;1 0 1;1 1 0]./norm([1 1 0]);%unit perfect vector matrix
numperf=size(perfmatrix,1); 
% calculate the power dissipated by the connected geometry and not splitting
        
for i=1:numtang
    glideplane(i,:)=cross(burgv(i,:),tangent(i,:));%glideplane of the perfect dislocation
    if(norm(glideplane(i,:))>eps)
        glideplane(i,:)=glideplane(i,:)/norm(glideplane(i,:));
    else
        glideplane(i,:)=links(linkid(i),4:6)/norm(links(linkid(i),4:6));
    end
    posnodeinsegment(i)=connectivity(nodeid,2*i+1);%This is the position of the node in the segment, can be 1 or 2
    othernodeid(i)=links(linkid(i),3-posnodeinsegment(i));
end
nodelist=[nodeid othernodeid]';
numnodes=length(nodelist);
cmax=size(connectivity,2);
clist=zeros(numnodes,cmax);
clist(1:numnodes,1)=[connectivity(nodelist,1)];
doesmove=1;
for k=1:numnodes
    clist(k,1)=connectivity(nodelist(k),1);
    clist(k,2:1+connectivity(nodelist(k),1))= linspace(1,connectivity(nodelist(k),1),connectivity(nodelist(k),1));
    if(rn(nodelist(k):flagpos)~=0)
        doesmove=0;
    end
end
[vntmp,fntmp]=feval(mobility,fseg,rn,links,connectivity,nodelist,clist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad); 
for k=1:numnodes
    Powermax=Powermax+dot(vntmp(k,:),fntmp(k,:));
end

Powermax=1.01*Powermax;%In Powertest(1) is the energy of the system unchanged

%Check the splitting in partials in the glide plane


for j=1:c,
    %posnodeinsegment(j)=connectivity(nodeid,2*j+1);
    %othernodeid(j)=links(linkid(j),3-posnodeinsegment(j));
    normplaneperf(j,:)=links(linkid(j),6:8);
    node1(j)=links(linkid(j),1);%vector containing the 1st nodes of each segment
    node2(j)=links(linkid(j),2);%vector containing the 2nd nodes of each segment 
end
newnodeid=size(rn,1)+1;% add new node
rn=[rn;zeros(1,flagpos)];
rn(newnodeid,:)=rn(nodeid,:);% with the same position as the node i
                                
                
for j=1:c,
    newlinkid(j)=size(links,1)+1;% in this case we have to add 2 new links
              
    links=[links;zeros(1,columlinks)];%add a new row to the link matrix
    stackforcevec=[stackforcevec;0 0 0]; 
                    
    if(posnodeinsegment(j)==1)
        links(newlinkid(j),1)=newnodeid;
        links(newlinkid(j),2)=node2(j);
    elseif(posnodeinsegment(j)==2)
        links(newlinkid(j),1)=node1(j);
        links(newlinkid(j),2)=newnodeid;
    end
                    
end
rnnewpos=rn;
linksnewpos=links;
stackforcevecnewpos=stackforcevec;
burgpart1=[];
burgpart2=[];
[burgpart1,burgpart2,burgpartcsp1,burgpartcsp2,splitdirgp,splitdircsp,crossslipplane]=burgervecpartialsnew(burgv,tangent,normplaneperf);%This function returns the burgers vectors of the partials depending on the perfect
if((norm(burgpart1)>eps)&(norm(burgpart2)>eps)&(norm(burgpartcsp1)>eps)&(norm(burgpartcsp2)>eps)&(norm(splitdirgp)>eps)&(norm(splitdircsp)>eps)&(norm(crossslipplane)>eps)&(doesmove==1))

    splitdirmatrix=[splitdirgp;-splitdirgp;splitdircsp;-splitdircsp];
    burg1matrix=[burgpart1;burgpart1;burgpartcsp1;burgpartcsp1];
    burg2matrix=[burgpart2;burgpart2;burgpartcsp2;burgpartcsp2];
    normalplanematrix=[normplaneperf;crossslipplane;crossslipplane];
    numsplitdir=size(splitdirmatrix,1);
    for i=1:numsplitdir
        
        Powertest=0;
        rn=rnnewpos;
        links=linksnewpos;
        stackforcevec=stackforcevecnewpos;
        noneedrd=0;
        rightcstan=0;
        burg1=burg1matrix(i,:);
        burg2=burg2matrix(i,:);
        glideplane=normalplanematrix(i,:);
        splitdir=splitdirmatrix(i,:);
        if((norm(burg1)~=0)&(norm(burg2)~=0)&(norm(glideplane)~=0))
            glideplane=glideplane/norm(glideplane);
%             if(burg1==burgpartcsp1)
%                 if(norm(cross(tangent(1,:)/norm(tangent(1,:)),tangent(2,:)/norm(tangent(2,:))))<rntol)
%                     for k=1:numperf
%                         if((norm(cross(tangent(1,:)/norm(tangent(1,:)),perfmatrix(k,:)))<rntol)&(norm(cross(tangent(2,:)/norm(tangent(2,:)),perfmatrix(k,:)))<rntol))
%                             if((posnodeinsegment(1)==2)&(dot(tangent(1,:),perfmatrix(k,:))>0))
%                                 rn(newnodeid,1:3)=rn(newnodeid,1:3)+(dot(tangent(1,:),perfmatrix(k,:))*perfmatrix(k,:))-tangent(1,:);
%                                 rn(nodeid,1:3)=rn(nodeid,1:3)+(dot(tangent(1,:),perfmatrix(k,:))*perfmatrix(k,:))-tangent(1,:);
%                             elseif((((posnodeinsegment(1)==1)&(dot(tangent(1,:),perfmatrix(k,:))<0)))|((posnodeinsegment(1)==2)&(dot(tangent(1,:),perfmatrix(k,:))<0)))
%                                 rn(newnodeid,1:3)=rn(newnodeid,1:3)+(dot(tangent(1,:),perfmatrix(k,:))*perfmatrix(k,:))+tangent(1,:);
%                                 rn(nodeid,1:3)=rn(nodeid,1:3)+(dot(tangent(1,:),perfmatrix(k,:))*perfmatrix(k,:))+tangent(1,:);
% %                             elseif((posnodeinsegment(1)==1)&(dot(tangent(1,:),perfmatrix(k,:))+0))
% %                                 rn(newnodeid,1:3)=rn(newnodeid,1:3)-(dot(tangent(1,:),perfmatrix(k,:))*perfmatrix(k,:))+tangent(1,:);
% %                                 rn(nodeid,1:3)=rn(nodeid,1:3)-(dot(tangent(1,:),perfmatrix(k,:))*perfmatrix(k,:))+tangent(1,:);
%                             end
%                             rightcstan=1;
%                             break;
%                         end
%                     end
%                 end
%             elseif(burg1==burgpart1)
%                 noneedrd=1;
%             end
%             if~((noneedrd==1)|(rightcstan==1))
%                 continue;
%             end
            for j=1:c,
                links(newlinkid(j),3:5)=burg1;%assign one of the Burgers vectors of the partials to the new link
                links(linkid(j),3:5)=burg2;%assign the other to the old link
                links(linkid(j),6:8)=glideplane;
                links(newlinkid(j),6:8)=glideplane;%assign the same normal as in the perfect case
            end
        
            rn(newnodeid,1:3)=rn(newnodeid,1:3)+((4*rann).*splitdir);%I choose this sense, but it doesn't matter
            rn(nodeid,1:3)=rn(nodeid,1:3)-((4*rann).*splitdir);
   
    
                % Now we have to set the direccion of the stacking fault force on the new partials
            for j=1:c,
                if(posnodeinsegment(j)==1)
                    tangpart1=(rn(othernodeid(j),1:3)-rn(newnodeid,1:3))/norm((rn(othernodeid(j),1:3)-rn(newnodeid,1:3)));%tangent unit to the new partial
                    L1=norm(rn(othernodeid(j),1:3)-rn(newnodeid,1:3));
                    normalnew1=cross(tangpart1,splitdir)/norm(cross(tangpart1,splitdir));%This is the flag we have to add to links so that this vector x tp gives the stacking fault force sense
                                              %The minus is because of the position of the node in the segment, because is the 1st
                                                                                
                    stackforcevec(newlinkid(j),:)=gamma.*normalnew1;%Assign the flag to the link
          %Now for the other partial
                    
                    tangpart2=(rn(othernodeid(j),1:3)-rn(nodeid,1:3))/norm((rn(othernodeid(j),1:3)-rn(nodeid,1:3)));%tangent unit to the new partial
                    L2=norm(rn(othernodeid(j),1:3)-rn(nodeid,1:3));
                    normalnew2=cross(tangpart2,-splitdir)/norm(cross(tangpart2,-splitdir));%This is the flag we have to add to links so that this vector x tp gives the stacking fault force sense
                                                       %The minus is because of the position of the node in the segment, because is the 1st
                        
                    stackforcevec(linkid(j),:)=gamma.*normalnew2;%Assign the flag to the link
                elseif(posnodeinsegment(j)==2)
                    tangpart1=-(rn(othernodeid(j),1:3)-rn(newnodeid,1:3))/norm((rn(othernodeid(j),1:3)-rn(newnodeid,1:3)));%tangent unit to the new partial
                    L1=norm(rn(othernodeid(j),1:3)-rn(newnodeid,1:3));
                    normalnew1=cross(tangpart1,splitdir)/norm(cross(tangpart1,splitdir));%This is the flag we have to add to links so that this vector x tp gives the stacking fault force sense
                                              %The minus is because of the position of the node in the segment, because is the 1st
                                                                                
                    stackforcevec(newlinkid(j),:)=gamma.*normalnew1;%Assign the flag to the link
          %Now for the other partial
                    
                    tangpart2=-(rn(othernodeid(j),1:3)-rn(nodeid,1:3))/norm((rn(othernodeid(j),1:3)-rn(nodeid,1:3)));%tangent unit to the new partial
                    L2=norm(rn(othernodeid(j),1:3)-rn(nodeid,1:3));
                    normalnew2=cross(tangpart2,-splitdir)/norm(cross(tangpart2,-splitdir));%This is the flag we have to add to links so that this vector x tp gives the stacking fault force sense
                                                       %The minus is because of the position of the node in the segment, because is the 1st
                        
                    stackforcevec(linkid(j),:)=gamma.*normalnew2;%Assign the flag to the link
                end
                
         %       if((i==3))&&(doSFT~=0)%
%                 if((i==3)||(i==4))&&(doSFT~=0)
%                    % Powertest=1e8;
%                     stackforcevec(linkid(j),:)=-stackforcevec(linkid(j),:);
%                     stackforcevec(newlinkid(j),:)=-stackforcevec(newlinkid(j),:);
%                 end 
            end
            
            
        
            [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
            linkslist=[linkid newlinkid]';
            fseg=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 flagpos]),links,appliedstress,linkslist,mobility,dopartials,stackforcevec,doinclusion,inclusion_pos_rad,doSFT,SFT_plane,doshielding);
            nodelist=[nodeid othernodeid newnodeid]';
            numnodes=length(nodelist);
            cmax=size(connectivity,2);
            clist=zeros(numnodes,cmax);
            clist(1:numnodes,1)=[connectivity(nodelist,1)];
            for k=1:numnodes
                clist(k,1)=connectivity(nodelist(k),1);
                clist(k,2:1+connectivity(nodelist(k),1))= linspace(1,connectivity(nodelist(k),1),connectivity(nodelist(k),1));
            end
            [vntmp,fntmp]=feval(mobility,fseg,rn,links,connectivity,nodelist,clist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad); 
            vec1=(rn(newnodeid,1:3)-rn(othernodeid(1),1:3))/norm(rn(newnodeid,1:3)-rn(othernodeid(1),1:3));
            vec2=(rn(newnodeid,1:3)-rn(othernodeid(2),1:3))/norm(rn(newnodeid,1:3)-rn(othernodeid(2),1:3));
            proj1=dot(vec1,rn(newnodeid,4:6));
            proj2=dot(vec2,rn(newnodeid,4:6));
            vec3=(rn(nodeid,1:3)-rn(othernodeid(1),1:3))/norm(rn(nodeid,1:3)-rn(othernodeid(1),1:3));
            vec4=(rn(nodeid,1:3)-rn(othernodeid(2),1:3))/norm(rn(nodeid,1:3)-rn(othernodeid(2),1:3));
            proj3=dot(rn(nodeid,4:6),vec3);
            proj4=dot(rn(nodeid,4:6),vec4);
            vec5=(rn(othernodeid(1),1:3)-rn(nodeid,1:3))/norm(rn(othernodeid(1),1:3)-rn(nodeid,1:3));
            vec6=(rn(othernodeid(1),1:3)-rn(newnodeid,1:3))/norm(rn(othernodeid(1),1:3)-rn(newnodeid,1:3));
            proj5=dot(rn(othernodeid(1),4:6),vec5);
            proj6=dot(rn(othernodeid(1),4:6),vec6);
            vec7=(rn(othernodeid(2),1:3)-rn(nodeid,1:3))/norm(rn(othernodeid(2),1:3)-rn(nodeid,1:3));
            vec8=(rn(othernodeid(2),1:3)-rn(newnodeid,1:3))/norm(rn(othernodeid(2),1:3)-rn(newnodeid,1:3));
            proj7=dot(rn(othernodeid(2),4:6),vec5);
            proj8=dot(rn(othernodeid(2),4:6),vec6);
            darea2dt=proj1+proj2+proj3+proj4+proj5+proj6;
            Powertest=0;       
            
            if(darea2dt>0)
            %if((proj1+proj2>0)&(proj3+proj4>0)&(proj5+proj6>0)&(proj7+proj8>0))
                
                for k=1:numnodes
                    Powertest=Powertest+dot(vntmp(k,:),fntmp(k,:));
                    rn(nodelist(k),4:6)=vntmp(k,:);
                end
            end
            if(Powertest>Powermax)
    
                rnnew=rn;
                linksnew=links;
                fsegnew=fseg;
                stackforcevecnew=stackforcevec;
                [connectivitynew,linksinconnectnew]=genconnectivity(rn,links,maxconnections);
                Powermax=Powertest;
                rn=rnnewpos;
                links=linksnewpos;
            %fseg=fsegold;
                stackforcevec=stackforcevecnewpos;
                reset=1;
                disp('splitting perfect');
            else
                rn=rnnewpos;
                links=linksnewpos;
            %fseg=fsegold;
                stackforcevec=stackforcevecnewpos;
            end
        end
    end
end
