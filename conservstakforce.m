function conservstakforce(stackforcevec,links,rn,connectivity)

[lrn,lrn2]=size(rn);
eps=1e-6;

for i=1:lrn
    c=connectivity(i,1);%number of connections of node i
    constackforce=[0 0 0];
    for j=1:c
        posnode=connectivity(i,2*j+1);
        if(posnode==1)
            constackforce=constackforce+stackforcevec(connectivity(i,2*j),:);
        elseif(posnode==2)
            constackforce=constackforce-stackforcevec(connectivity(i,2*j),:);
        end
    end
    if(norm(constackforce)>eps)&(rn(i,lrn2)==0)
       
        disp(sprintf('FAIL IN STACKING FAULT FORCE VECTOR CONSERVATION IN NODE %d',i));
        pause
    end
end