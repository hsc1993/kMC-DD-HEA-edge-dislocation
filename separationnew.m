function [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,stackforcevecnew]=separationnew(rn,links,connectivity,linksinconnect,fseg,mobility,MU,NU,a,Ec,mindist,appliedstress,dopartials,stackforcevec,rann,gamma,rntol,maxconnections,doSFT,SFT_plane,doinclusion,inclusion_pos_rad,doshielding,SurfacePlane,dobox,boxD)
[lrn,lrn2]=size(rn);
eps=1e-6;
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
normstairrodobtuseSF=norm(gamma.*(1/sqrt(3))*[2 0 0]);
normHirthacuteSF=norm(gamma.*(1/sqrt(3))*[2 2 0]);
normstairrodBV=norm((1/6).*[1 1 0]);
normHirthBV=norm((1/3).*[1 0 0]);

% search through connectivity list
for i=1:lrn
    c=connectivity(i,1);
    if (c>3) & (rn(i,lrn2)==0)
        % disp(sprintf('separation: node %d has %d arms',i,c));
        % a high order node has been found
        % calculate the power dissipated by the connected geometry and not splitting
        ft=zeros(1,3);
        for j=1:c
            linkid(j)=connectivity(i,2*j);
            pos(j)=connectivity(i,2*j+1);
            othernodeid(j)=links(linkid(j),3-pos(j));
            
        end
        nodelist=[i othernodeid]';
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
        for k=1:numnodes
            Powermax=Powermax+dot(vntmp(k,:),fntmp(k,:));
        end
        Powermax=1.01*Powermax;
        %initialize the splitting mode that will be undertaken
        splittingmode=0;
        
        % build the list of possible splittingmodes of the multinode
        conlist=buildsplitmodelist(c);
        numsplitmodes=size(conlist,1);
        % conlist is a matrix with dimension numsplitting modes by number of connections
        % it consists of a series of ones and twos detailing which connection in the original
        % connectivity list is connected in the temporary dislocation structure that is being considered
        
        % save the current configuration of the node to be considered for splitting
        refposveli=rn(i,1:6);
        refconnecti=connectivity(i,:);
        for j=1:c
            reffsegi(j,:)=fseg(refconnecti(2*j),:);
        end
        reflinksinconnect=linksinconnect;
        refrn=rn;
        refconnectivity=connectivity;
        reffseg=fseg;
        reflinks=links;
        refstackforcevec=stackforcevec;
        % begin investigating the power dissipation of the separated configurations
        for j=1:numsplitmodes
            s=conlist(j,:);
            for k=1:c
                count=c-k+1;
                if s(count)==1
                    s(count)=count;
                else
                    s(count)=[];
                end
            end
            [rn,links,connectivity,linksinconnect,stackforcevec]=splitnode(rn,links,connectivity,linksinconnect,i,s,refposveli,stackforcevec,0);
            %lastlink=size(links,1);
            %if(((abs(norm(links(lastlink,3:5))-normstairrodBV)<eps)&(abs(norm(stackforcevec(lastlink,:))-normstairrodobtuseSF)<eps))|((abs(norm(links(lastlink,3:5))-normHirthBV)<eps)&(abs(norm(stackforcevec(lastlink,:))-normHirthacuteSF)<eps)))
            %    rn=refrn;
            %    links=reflinks;
            %    stackforcevec=refstackforcevec;
            %    fseg=reffseg;
            %    connectivity=refconnectivity;
            %    linksinconnect=reflinksinconnect;
            %    continue;
            %end
            lastnode=size(rn,1);
            nodelist=[i lastnode othernodeid]';
            numnodes=length(nodelist);
            cmax=size(connectivity,2);
            clist=zeros(numnodes,cmax);
            clist(1:numnodes,1)=[connectivity(nodelist,1)];
            for k=1:numnodes
                  %clist(k,1)=connectivity(nodelist(k),1);
                clist(k,2:1+connectivity(nodelist(k),1))= linspace(1,connectivity(nodelist(k),1),connectivity(nodelist(k),1));
            end 
             
            [vntmp,fntmp]=feval(mobility,fseg,rn,links,connectivity,nodelist,clist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad);
            
            
            vd1=vntmp(1,:)*vntmp(1,:)';
            vd2=vntmp(2,:)*vntmp(2,:)';
         %   if(dot(vntmp(1,:),vntmp(2,:))>0)
            
                if vd1>=vd2
                    dir=-vntmp(1,:)./sqrt(vd1);
                    rn(i,1:3)=rn(i,1:3)-(2*mindist)*dir;
                else
                    dir=vntmp(2,:)./sqrt(vd2);
                    rn(lastnode,1:3)=rn(lastnode,1:3)+(2*mindist)*dir;
                end
                linkslen=size(links,1);
                if(norm(links(linkslen,6:8))<eps)
                    burgnew=links(linkslen,3:5);
                    links(linkslen,6:8)=cross(dir,burgnew)/norm(cross(dir,burgnew));
                end
                fseg=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,appliedstress,[],mobility,dopartials,stackforcevec,doinclusion,inclusion_pos_rad,doSFT,SFT_plane,doshielding,SurfacePlane,dobox,boxD);   
                [vntmp,fntmp]=feval(mobility,fseg,rn,links,connectivity,nodelist,clist,mobility,dopartials,stackforcevec,rntol,doSFT,SFT_plane,doinclusion,inclusion_pos_rad);
                rn(nodelist,4:6)=vntmp;
                vdiff=vntmp(2,:)-vntmp(1,:);
                if vdiff*dir'>0
                    Powertest=0;
                    for k=1:numnodes
                        Powertest=Powertest+fntmp(k,:)*vntmp(k,:)';
                    end
                    if Powermax<Powertest
                        Powermax=Powertest;
                        rnnew=rn;
                        linksnew=links;
                        stackforcevecnew=stackforcevec;
                        fsegnew=fseg;
                        [connectivitynew,linksinconnectnew]=genconnectivity(rn,links,maxconnections);
                        Powermax=Powertest;
                        
                        disp('splitting is taken place');
                                           
                    end
                end
                
                %  end
            rn=rnold;
            links=linksold;
            stackforcevec=stackforcevecold;
            fseg=fsegold;
            [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
        end
    end
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [conlist]=buildsplitmodelist(numconnections)
    c=numconnections;
    maxsplitarms=floor(0.5*c);
    % calculate the number of possible splitting modes
    numsplitmodes=0;
    conlist=[];
    count=1;
    for j=2:maxsplitarms
        subnumber=factorial(c)/factorial(c-j)/factorial(j);
        splitmodes=createsplitmodelist(c,j);
        if (2*j==c)% this is a special case that handles mirror symmetry situations
            subnumber=0.5*subnumber;
        end
        for k=1:subnumber
            conlist(count,:)=zeros(1,c);
            for i=1:j
                conlist(count,splitmodes(k,i))=1;
            end
            count=count+1;
        end
        numsplitmodes=numsplitmodes+subnumber;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

function [splitmodes]=createsplitmodelist(nodeconnectivity,numsplitarms)
%initialize the size of splitmodes and fill with zeros
    splitmodes=zeros( factorial(nodeconnectivity) / factorial(nodeconnectivity-numsplitarms) / factorial(numsplitarms) , numsplitarms );
    if (numsplitarms==1)
        % if numsplitarms is 1 then the answer is simply a count from one to nodeconnectivity
        splitmodes=linspace(1,nodeconnectivity,nodeconnectivity)';
    else
        % if numsplitarms is greater than one go into a recursive algorithm to build up the system
        lastone=0;
        for i=1:nodeconnectivity-numsplitarms+1
            % initialize the subslitmodes matrix b
            addlength= factorial(nodeconnectivity-i) / factorial(nodeconnectivity-i-numsplitarms+1) / factorial(numsplitarms-1);
            b=zeros(addlength,numsplitarms-1);
            % calculate the subsplitmodes with a call to splitmodes
            b=createsplitmodelist(nodeconnectivity-i,numsplitarms-1);
            % add the subspitmode list b to the growing splitmode list
            splitmodes(lastone+1:lastone+addlength,1)=i.*ones(addlength,1);
            splitmodes(lastone+1:lastone+addlength,2:numsplitarms)=i+b( 1:addlength , 1:numsplitarms-1 );
            % update the current position of the last valid line of the splitmode system
            lastone=lastone+addlength;
        end       
    end


%total number of splitting options
%4 arm multinode
%may split into 3 two arm splitting modes      0.5 * n! / ( (n-2)! * 2! )       3 total
%5 arm mulitnode
%may split into 10 two arm splitting modes           n! / ( (n-2)! * 2! )      10 total
%6 arm mulitinode
%may split into 15 two arm splitting modes or        n! / ( (n-2)! * 2! )
%may split into 10 three arm spliting modes    0.5 * n! / ( (n-3)! * 3! )      25 total
%7 arm multinode
%may split into 21 two arm splitting modes or        n! / ( (n-2)! * 2! )
%may split into 35 thee arm splitting modes          n! / ( (n-3)! * 3! )      46 total
%8 arm multinode
%may split into 28 two arm splitting modes or        n! / ( (n-2)! * 2! )
%may split into 56 three arm splitting modes or      n! / ( (n-3)! * 3! )
%may split into 35 four arm splitting modes    0.5 * n! / ( (n-4)! * 4! )     119 total 
%9 arm multinode
%may split into 36 two arm splitting modes or        n! / ( (n-2)! * 2! )
%may split into 84 three arm splitting modes or      n! / ( (n-3)! * 3! )
%may split into 126 four arm splitting modes         n! / ( (n-4)! * 4! )     246 total
%10 arm multinode
%may split into 45 two arm splitting modes or        n! / ( (n-2)! * 2! )
%may split into 120 three arm splitting modes or     n! / ( (n-3)! * 3! )
%may split into 210 four arm splitting modes or      n! / ( (n-4)! * 4! )
%may split into 126 five arm splitting modes   0.5 * n! / ( (n-5)! * 5! )     501 total
% the division by 2 for some of the cases is done outside this function call when appropriate.    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


