%This function handles the general multinode splitting

function [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew,printnode]=crossmultinode(rn,links,connectivity,linksinconnect,fseg,stackforcevec,burgv,tangent,nodeid,linkid,mobility,dopartials,rann,gamma,MU,NU,a,Ec,appliedstress,rntol,areamin,lmin,printnode,doSFT,SFT_plane,doinclusion,inclusion_pos_rad)

%rntol=0.1;
numarms=size(tangent,1);
eps=1e-3;
printnode=printnode;
crosstol=0.1;%areamin;%THERE ARE SOME PROBLEMS WITH THIS VALUE. CHECKKKKKKKKKKK
maxconnections=8;%Can take it from the input
delta=1e-16;
perfmatrix=[1 -1 0;-1 0 1;0 -1 1;0 1 1;1 0 1;1 1 0]./norm([1 1 0]);
numperfdir=size(perfmatrix,1);
normperfburg=norm(1/2.*[1 1 0]);
normpartburg=norm(1/6.*[1 1 2]);
normstairrod=norm(1/6.*[1 1 0]);
crossplane=[1 1 1; -1 1 -1; -1 -1 1; 1 -1 -1]/norm([1 1 1]);
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
[numnodes,flagpos]=size(rn);%Total number of nodes
[numlinks,columlinks]=size(links);%Total number of links
[numfseg,columfseg]=size(fseg);
[numstack,columstack]=size(stackforcevec);
for i=1:numarms
    posnodeinsegment(i)=connectivity(nodeid,2*i+1);%This is the position of the node in the segment, can be 1 or 2
    if (posnodeinsegment(i)==2)
        tangent(i,:)=-tangent(i,:);
    end
    othernodeid(i)=links(linkid(i),3-posnodeinsegment(i));
end
nodelist=[nodeid othernodeid]';
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
Powermax=1.01*Powermax;
%Powermax=0;%IF WE DON'T CHECK ALL THE NODES
for i=1:numarms
    iglideplane=links(linkid(i),6:8)/norm(links(linkid(i),6:8));
    for j=1:numarms
        if((i~=j)&(norm(cross(tangent(i,:),tangent(j,:)))>crosstol))
            jglideplane=links(linkid(j),6:8)/norm(links(linkid(j),6:8));
            newplane=cross(tangent(i,:),tangent(j,:))/norm(cross(tangent(i,:),tangent(j,:)));
           
            if ((norm(newplane-iglideplane)<eps)&(norm(newplane-jglideplane)<eps))%That means that both segments are in the same plane and that this plane is the glide plane of both
                if(abs(norm(burgv(i,:))-normperfburg)<eps)  %The only thing can happen is that a perfect dislocation splits into partials
                %[burgpart1,burgpart2]=burgervecpartials(burgv(i,:),tangent(i,:),iglideplane);
                %Pa luego
                end
            end
            numcrossplanes=size(crossplane,1);
            for m=1:numcrossplanes
                if((norm(cross(newplane,crossplane(m,:)))<crosstol)&(norm(iglideplane-newplane)>eps)) %If both segments form a [1 1 1] plane, and the glideplane of the segment changes
                    for k=1:numperfdir
                        for l=1:numperfdir
                            if((norm(cross(tangent(i,:)/norm(tangent(i,:)),perfmatrix(k,:)))<crosstol)&(norm(cross(tangent(j,:)/norm(tangent(j,:)),perfmatrix(l,:)))<crosstol)&(k~=l)&(norm(tangent(i,:))>lmin)&(norm(tangent(j,:))>lmin)) %The segments have to be in [1 1 0] direction. CHECKKKKKK
                                %if(abs(norm(burgv(i,:))-normperfburg)<eps)&(abs(dot(burgv(i,:)/normperfburg,tangent(i,:)/norm(tangent(i,:))))>0.2)  %The only thing can happen is that a perfect dislocation splits into partials
                                %[burgpart1,burgpart2,splitdir,crossplane]=burgervecpartials(burgv(i,:),tangent(i,:),iglideplane);
                                %[rn,links,stackforcevec,newnodeidi,newlinkidi]=addnewnode(rn,links,stackforcevec,tangent(i,:),linkid(i),nodeid,posnodeinsegment(i),rann,[]);
                                %[rn,links,stackforcevec,newnodeidj,newlinkidj]=addnewnode(rn,links,stackforcevec,tangent(j,:),linkid(j),nodeid,posnodeinsegment(j),rann,[]);
                                %newlinkadd=size(links,1)+1;
                                %burgpart1(2,:)=zeros(1,3);
                                %burgpart1(2,:)=burgpart2;
                                %burgpart2(2,:)=zeros(1,3);
                                %burgpart2(2,:)=burgpart1;
                                %possiblesplit=size(burgpart1,1);
                                %for p=1:possiblesplit
                                    %links=[links;zeros(1,columlinks)];
                                    %stackforcevec=[stackforcevec;zeros(1,3)];
                                    %if(posnodeinsegment(i)==1)
                                     %   links(newlinkadd,1)=newnodeidj;
                                     %   links(newlinkadd,2)=newnodeidi;
                                     %elseif(posnodeinsegment(i)==2)
                                     %   links(newlinkadd,1)=newnodeidi;
                                     %   links(newlinkadd,2)=newnodeidj;
                                     %end
                            
                                    %if(posnodeinsegment(i)==posnodeinsegment(j))
                                    %    newburgvec=newburpar+burgv(j,:);
                                    %else
                                    %    newburgvec=burgv(j,:)-newburpar;
                                    %end 
                                    %if((abs(norm(newburgvec)-normperfburg)<eps)|(abs(norm(newburgvec)-normpartburg)<eps)|(abs(norm(newburgvec)-normstairrod)<eps))
                                  
                                
                                    %end
                                if(abs(norm(burgv(i,:))-normpartburg)<eps)&(abs(dot(burgv(i,:)/normpartburg,tangent(i,:)/norm(tangent(i,:))))>0.2) %If it is a Shockley partial and it is not an edge
                                    [newburpar,newstairod,splitdir]=tablepartialstairod(burgv(i,:),tangent(i,:),rntol);
                                    if(norm(newburpar)>0)&(norm(newstairod)>0)
                                        [rn,links,stackforcevec,newnodeidi,newlinkidi]=addnewnode(rn,links,stackforcevec,tangent(i,:),linkid(i),nodeid,posnodeinsegment(i),rann,[]);
                                        [rn,links,stackforcevec,newnodeidj,newlinkidj]=addnewnode(rn,links,stackforcevec,tangent(j,:),linkid(j),nodeid,posnodeinsegment(j),rann,[]);
                                        newlinkadd=size(links,1)+1;
                                    
                                        links=[links;zeros(1,columlinks)];
                                        stackforcevec=[stackforcevec;zeros(1,3)];
                                        if(posnodeinsegment(i)==1)
                                            links(newlinkadd,1)=newnodeidj;
                                            links(newlinkadd,2)=newnodeidi;
                                        elseif(posnodeinsegment(i)==2)
                                            links(newlinkadd,1)=newnodeidi;
                                            links(newlinkadd,2)=newnodeidj;
                                        end
                            
                                        if(posnodeinsegment(i)==posnodeinsegment(j))
                                            newburgvec=burgv(j,:)+newburpar;
                                        else
                                            newburgvec=burgv(j,:)-newburpar;
                                        end 
                                    
                                        if((abs(norm(newburgvec)-normperfburg)<eps)|(abs(norm(newburgvec)-normpartburg)<eps)|(abs(norm(newburgvec)-normstairrod)<eps))
                               
                                            t1=(rn(links(newlinkadd,2),1:3)-rn(links(newlinkadd,1),1:3))/norm((rn(links(newlinkadd,2),1:3)-rn(links(newlinkadd,1),1:3)));
                                            t2=(rn(links(newlinkidi,2),1:3)-rn(links(newlinkidi,1),1:3))/norm((rn(newnodeidi,1:3)-rn(nodeid,1:3)));%Aux vector to give me the stacking fault direction
                                            t3=(rn(links(newlinkidj,2),1:3)-rn(links(newlinkidj,1),1:3))/norm((rn(newnodeidj,1:3)-rn(nodeid,1:3)));
                                            if(posnodeinsegment(i)==1)
                                                stackforcevec(newlinkadd,:)=gamma.*cross(t1,t2)/norm(cross(t1,t2));
                                                stackforcevec(newlinkidi,:)=stackforcevec(linkid(i),:)-stackforcevec(newlinkadd,:);
                                            else
                                                stackforcevec(newlinkadd,:)=gamma.*cross(t2,t1)/norm(cross(t2,t1)); 
                                                stackforcevec(newlinkidi,:)=stackforcevec(linkid(i),:)-stackforcevec(newlinkadd,:);
                                            end
                                            if(posnodeinsegment(i)==posnodeinsegment(j))
                                                stackforcevec(newlinkidj,:)=stackforcevec(newlinkadd,:)+stackforcevec(linkid(j),:);
                                            else
                                                stackforcevec(newlinkidj,:)=stackforcevec(linkid(j),:)-stackforcevec(newlinkadd,:);
                                            end 
                                            links(newlinkadd,3:5)=newburpar;
                                            links(newlinkadd,6:8)=newplane;
                                            links(newlinkidi,3:5)=newstairod;
                                            if(norm(stackforcevec(newlinkidi,:))>eps)
                                                links(newlinkidi,6:8)=stackforcevec(newlinkidi,:)/norm(stackforcevec(newlinkidi,:));
                                            else
                                                links(newlinkidi,6:8)=cross(links(newlinkidi,3:5),t2)/norm(cross(links(newlinkidi,3:5),t2));
                                            end
                                            links(newlinkidj,3:5)=newburgvec;
                                            if(norm(stackforcevec(newlinkidj,:))>eps)
                                                links(newlinkidj,6:8)=stackforcevec(newlinkidj,:)/norm(stackforcevec(newlinkidj,:));
                                            else
                                                links(newlinkidj,6:8)=cross(links(newlinkidj,3:5),t3)/norm(cross(links(newlinkidj,3:5),t3));
                                            end
                                            
                                            [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
                                            fseg=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 flagpos]),links,appliedstress,[],mobility,dopartials,stackforcevec,doinclusion,inclusion_pos_rad);    
                                            nodelist=[nodeid newnodeidi newnodeidj othernodeid]';
                                            numnodes=length(nodelist);
                                            cmax=size(connectivity,2);
                                            clist=zeros(numnodes,cmax);
                                            clist(1:numnodes,1)=[connectivity(nodelist,1)];
                                            for m=1:numnodes
                                            %clist(k,1)=connectivity(nodelist(k),1);
                                                clist(m,2:1+connectivity(nodelist(m),1))= linspace(1,connectivity(nodelist(m),1),connectivity(nodelist(m),1));
                                            end   
                                            [vntmp,fntmp]=feval(mobility,fseg,rn,links,connectivity,nodelist,clist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad);
                                            for m=1:numnodes
                                                rn(nodelist(m),4:6)=vntmp(m,:);
                                            end
                            %Check the area change
                    
                                            vec1=rn(newnodeidi,1:3)-rn(nodeid,1:3);
                                            vec2=rn(newnodeidj,1:3)-rn(nodeid,1:3);
                                            vec3=vec2-vec1;
                                            r1=sqrt(vec1*vec1');
                                            r2=sqrt(vec2*vec2');
                                            r3=sqrt(vec3*vec3');
                                            s=0.5*(r1+r2+r3);
                                            area2=(s*(s-r1)*(s-r2)*(s-r3));
                                            dvec1dt=rn(newnodeidi,4:6)-rn(nodeid,4:6);
                                            dvec2dt=rn(newnodeidj,4:6)-rn(nodeid,4:6);
                                            dvec3dt=dvec2dt-dvec1dt;
                                            dr1dt=(vec1*dvec1dt')/(r1+delta);
                                            dr2dt=(vec2*dvec2dt')/(r2+delta);
                                            dr3dt=(vec3*dvec3dt')/(r3+delta);
                                            dsdt=0.5*(dr1dt+dr2dt+dr3dt);
                                            darea2dt=dsdt*(s-r1)*(s-r2)*(s-r3);
                                            darea2dt=darea2dt+s*(dsdt-dr1dt)*(s-r2)*(s-r3);
                                            darea2dt=darea2dt+s*(s-r1)*(dsdt-dr2dt)*(s-r3);
                                            darea2dt=darea2dt+s*(s-r1)*(s-r2)*(dsdt-dr3dt);
                                            Powertest=0;
     
                                            if(darea2dt>0)
                                            %vnglide=zeros(numnodes,3);
                                                for m=1:numnodes
                                                    Powertest=Powertest+dot(vntmp(m,:),fntmp(m,:));
                                                end 
                                            end
                                            if(Powertest>Powermax)
                                                rnnew=rn;
                                                linksnew=links;
                                                stackforcevecnew=stackforcevec;
                                                fsegnew=fseg;
                                                [connectivitynew,linksinconnectnew]=genconnectivity(rn,links,maxconnections);
                                                Powermax=Powertest;
                                                rn=rnold;
                                                links=linksold;
                                                stackforcevec=stackforcevecold;
                                                fseg=fsegold;
                                                disp('multicross-slip is taken place');
                                                printnode=[printnode newnodeidi newnodeidj];
                                            end
                                        end
                                        rn=rnold;
                                        links=linksold;
                                        stackforcevec=stackforcevecold;
                                        fseg=fsegold;
                                        [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
                                    end
                                    
                                elseif(abs(norm(burgv(i,:))-normstairrod)<eps)%&(abs(dot(burgv(i,:)/normstairrod,tangent(i,:)/norm(tangent(i,:))))<rntol)%If it is a Stair-rod and it is an edge      
                                    [burgpart1,burgpart2,stacknewlink,stacklink,otherplane]=tablesplitstairrod(burgv(i,:),newplane,stackforcevec(linkid(i),:),crosstol,gamma);%Burgpart1 is gonna be the burgers vector of the newlinkadd 
                                    
                                    possiblesplit=size(burgpart1,1);
                                    for p=1:possiblesplit
                                        if(norm(burgpart1(p,:))>0)&(norm(burgpart2(p,:))>0)
                                            [rn,links,stackforcevec,newnodeidi,newlinkidi]=addnewnode(rn,links,stackforcevec,tangent(i,:),linkid(i),nodeid,posnodeinsegment(i),rann,[]);
                                            [rn,links,stackforcevec,newnodeidj,newlinkidj]=addnewnode(rn,links,stackforcevec,tangent(j,:),linkid(j),nodeid,posnodeinsegment(j),rann,[]);
                                            newlinkadd=size(links,1)+1;
                                            links=[links;zeros(1,columlinks)];
                                            stackforcevec=[stackforcevec;zeros(1,3)];
                                            if(posnodeinsegment(i)==1)
                                                links(newlinkadd,1)=newnodeidj;
                                                links(newlinkadd,2)=newnodeidi;
                                            elseif(posnodeinsegment(i)==2)
                                                links(newlinkadd,1)=newnodeidi;
                                                links(newlinkadd,2)=newnodeidj;
                                            end
                            
                                            if(posnodeinsegment(i)==posnodeinsegment(j))
                                                newburgvec=burgpart1(p,:)+burgv(j,:);
                                            else
                                                newburgvec=burgv(j,:)-burgpart1(p,:);
                                            end 
                                            if((abs(norm(newburgvec)-normperfburg)<eps)|(abs(norm(newburgvec)-normpartburg)<eps)|(abs(norm(newburgvec)-normstairrod)<eps))
                                                stackforcevec(newlinkadd,:)=stacknewlink;%dot(stackforcevec(linkid(i),:),newplane).*newplane;
                                                stackforcevec(newlinkidi,:)=stacklink;%stackforcevec(linkid(i),:)-stackforcevec(newlinkadd,:);

                            
                                                if(posnodeinsegment(i)==posnodeinsegment(j))
                                                    stackforcevec(newlinkidj,:)=stackforcevec(newlinkadd,:)+stackforcevec(linkid(j),:);
                                                else
                                                    stackforcevec(newlinkidj,:)=stackforcevec(linkid(j),:)-stackforcevec(newlinkadd,:);
                                                end 
                                                links(newlinkadd,3:5)=burgpart1(p,:);
                                                links(newlinkadd,6:8)=newplane;
                                                links(newlinkidi,3:5)=burgpart2(p,:);
                                                links(newlinkidi,6:8)=otherplane;%cross(tangent(i,:),burgpart2(p,:))/norm(cross(tangent(i,:),burgpart2(p,:)));
                                                links(newlinkidj,3:5)=newburgvec; 
                                                if(norm(cross(tangent(j,:),newburgvec))>eps)
                                                    links(newlinkidj,6:8)=cross(tangent(j,:),newburgvec)/norm(cross(tangent(j,:),newburgvec));
                                                elseif(norm(stackforcevec(newlinkidj,:))>eps)
                                                    links(newlinkidj,6:8)=stackforcevec(newlinkidj,:)/norm(stackforcevec(newlinkidj,:));
                                                else %The segment is screw perfect
                                                    [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
                                                    fseg=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 flagpos]),links,appliedstress,[],mobility,dopartials,stackforcevec);
                                                    cmax=size(connectivity,2);
                                                    clist=zeros(1,cmax);
                                                    clist(1,1)=[connectivity(newnodeidj,1)];
                                                    clist(1,2:1+connectivity(newnodeidj,1))= linspace(1,connectivity(newnodeidj,1),connectivity(newnodeidj,1));
                                                    [vntmp,fntmp]=feval(mobility,fseg,rn,links,connectivity,newnodeidj,clist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane);
                                                    links(newlinkidj,6:8)=cross(tangent(j,:),vntmp)/norm(cross(tangent(j,:),vntmp));
                                                end
                                                [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
                                                fseg=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 flagpos]),links,appliedstress,[],mobility,dopartials,stackforcevec,doinclusion,inclusion_pos_rad);    
                                                nodelist=[nodeid newnodeidi newnodeidj othernodeid]';
                                                numnodes=length(nodelist);
                                                cmax=size(connectivity,2);
                                                clist=zeros(numnodes,cmax);
                                                clist(1:numnodes,1)=[connectivity(nodelist,1)];
                                                for m=1:numnodes
                                                %clist(k,1)=connectivity(nodelist(k),1);
                                                    clist(m,2:1+connectivity(nodelist(m),1))= linspace(1,connectivity(nodelist(m),1),connectivity(nodelist(m),1));
                                                end   
                                                [vntmp,fntmp]=feval(mobility,fseg,rn,links,connectivity,nodelist,clist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad);
                                                rn(nodelist,4:6)=vntmp;
                                                                        %Check the area change
                    
                                                vec1=rn(newnodeidj,1:3)-rn(nodeid,1:3);
                                                vec2=rn(newnodeidi,1:3)-rn(nodeid,1:3);
                                                vec3=vec2-vec1;
                                                r1=sqrt(vec1*vec1');
                                                r2=sqrt(vec2*vec2');
                                                r3=sqrt(vec3*vec3');
                                                s=0.5*(r1+r2+r3);
                                                area2=(s*(s-r1)*(s-r2)*(s-r3));
                                                dvec1dt=rn(newnodeidj,4:6)-rn(nodeid,4:6);
                                                dvec2dt=rn(newnodeidi,4:6)-rn(nodeid,4:6);
                                                dvec3dt=dvec2dt-dvec1dt;
                                                dr1dt=(vec1*dvec1dt')/(r1+delta);
                                                dr2dt=(vec2*dvec2dt')/(r2+delta);
                                                dr3dt=(vec3*dvec3dt')/(r3+delta);
                                                dsdt=0.5*(dr1dt+dr2dt+dr3dt);
                                                darea2dt=dsdt*(s-r1)*(s-r2)*(s-r3);
                                                darea2dt=darea2dt+s*(dsdt-dr1dt)*(s-r2)*(s-r3);
                                                darea2dt=darea2dt+s*(s-r1)*(dsdt-dr2dt)*(s-r3);
                                                darea2dt=darea2dt+s*(s-r1)*(s-r2)*(dsdt-dr3dt);
                                                Powertest=0;
     
                                                if(darea2dt>0)
                                               % vnglide=zeros(numnodes,3);
                                                    for m=1:numnodes
                                                        Powertest=Powertest+dot(vntmp(m,:),fntmp(m,:));
                                                    end 
                                                end
                                                if(Powertest>Powermax)
                                                    rnnew=rn;
                                                    linksnew=links;
                                                    stackforcevecnew=stackforcevec;
                                                    fsegnew=fseg;
                                                    [connectivitynew,linksinconnectnew]=genconnectivity(rn,links,maxconnections);
                                                    Powermax=Powertest;
                                                    rn=rnold;
                                                    links=linksold;
                                                    stackforcevec=stackforcevecold;
                                                    fseg=fsegold;
                                                    disp('multicross-slip is taken place');
                                                    printnode=[printnode newnodeidi newnodeidj];
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
                    end
                end
            end
        end
    end
end
                                
                                
                                
                            
                
    