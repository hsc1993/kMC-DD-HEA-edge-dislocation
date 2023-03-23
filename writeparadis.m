function writeparadis(fname,rn,links,MU,NU,a)
%write rn, links into ParaDiS restart file format
%rn has to be a compact array (no deleted nodes in the middle)

%settings
BOXL=500;
minSideX = -BOXL;
maxSideX =  BOXL;
minSideY = -BOXL;
maxSideY =  BOXL;
minSideZ = -BOXL;
maxSideZ =  BOXL;
xBoundType = 1; %Free
yBoundType = 1; %Free
zBoundType = 1; %Free
maxSeg = 4000.;  %to disable remesh
minSeg =  0.01;  %to disable remesh
readRijmfile = 0; %no PBC correction
subdiv=4;

%internal data structure
[NMAX,m]=size(rn);
[LINKMAX,m]=size(links);

%number of nodes
np=sum((rn(:,4)~=-1));

%build link list
list=zeros(np,100);
for j=1:LINKMAX,
    n0=links(j,1);
    n1=links(j,2);
    if(n0~=0)&(n1~=0)
        list(n0,1)=list(n0,1)+1;
        list(n0,list(n0,1)+1)=j;
        list(n1,1)=list(n1,1)+1;
        list(n1,list(n1,1)+1)=j;
    end
end


%write text restart file
fp=fopen(fname,'w');

fprintf(fp,'###ParaDiS restart file created by dd3d.m (use -*-shell-script-*- format)\n\n');
fprintf(fp,'rc = %f\n\n',a);
fprintf(fp,'subdiv = %d\n',subdiv);
fprintf(fp,'shearModulus = %f\n',MU);
fprintf(fp,'pois = %f\n',NU);
fprintf(fp,'minSideX = %f\n',minSideX);
fprintf(fp,'maxSideX = %f\n',maxSideX);
fprintf(fp,'minSideY = %f\n',minSideY);
fprintf(fp,'maxSideY = %f\n',maxSideY);
fprintf(fp,'minSideZ = %f\n',minSideZ);
fprintf(fp,'maxSideZ = %f\n',maxSideZ);

fprintf(fp,'xBoundType = %d\n',xBoundType);
fprintf(fp,'yBoundType = %d\n',yBoundType);
fprintf(fp,'zBoundType = %d\n',zBoundType);

fprintf(fp,'maxSeg = %f\n',maxSeg);
fprintf(fp,'minSeg = %f\n',minSeg);

fprintf(fp,'readRijmfile = %d\n',readRijmfile);

fprintf(fp,'\nconfig = [\n');
fprintf(fp,'%f %f %f\n',minSideX,minSideY,minSideZ);
fprintf(fp,'%f %f %f\n',maxSideX,maxSideY,maxSideZ);
fprintf(fp,'#Burgers vector array (0: disabled)\n0\n');    
fprintf(fp,'#number of nodes\n%d\n',np);
fprintf(fp,'### (Primary lines: node_id, old_id, x, y, z, numNbrs, constraint, domain, index)\n');
fprintf(fp,'### (Secondary lines:   nbr[i], burgX[i] burgY[i] burgZ[i] nx[i] ny[i] nz[i])\n');

weight=0; domainID=0;
for i=1:np,
    numNbrs=list(i,1);
    fprintf(fp, '  %7d %6d %14.4f %14.4f %14.4f %3d %3d %4d %5d\n',...
        i,weight,rn(i,1),rn(i,2),rn(i,3),numNbrs,rn(i,4),domainID,i-1);
    for j=1:numNbrs,
        k=list(i,j+1);
        n0=links(k,1);
        n1=links(k,2);
        bv=links(k,3:5)';
        nv=links(k,6:8)';
        if(n0==i)
            nb=n1;
        else
            nb=n0;
            bv=-bv;
        end
        fprintf(fp, '        %6d %16.10f %16.10f %16.10f\n               %16.10f %16.10f %16.10f\n',...
            nb,bv(1),bv(2),bv(3),nv(1),nv(2),nv(3));
    end
end

fprintf(fp,'###domain boundaries (nXdoms, nYdoms, nZdoms)\n');
fprintf(fp,'   1    1    1\n');
fprintf(fp,'###domain boundaries (X,   Y,    Z)\n');
fprintf(fp,'   %f\n',minSideX);
fprintf(fp,'           %f\n',minSideY);
fprintf(fp,'                   %f\n',minSideZ);
fprintf(fp,'                   %f\n',maxSideZ);
fprintf(fp,'           %f\n',maxSideY);
fprintf(fp,'   %f\n',maxSideX);

fprintf(fp,']\n');

fclose(fp);
