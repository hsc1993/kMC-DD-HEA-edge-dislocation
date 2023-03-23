%This function gives me the burgers vectors in case a stair-rod splits

function [burg1,burg2,splitdir,otherplane]=tablefrankpartialsnew(burgv,tangent,rntol,stackforcevec)

burg1=zeros(5,3);
burg2=zeros(5,3);
otherplane=zeros(5,3);
splitdir=zeros(5,3);
rntol=1e-3;

if(norm(burgv-(1/3.*[-1 1 -1]))<rntol)  %aA frank partial
    burg1(3,:)=1/6.*[-2 -1 1];
    burg2(3,:)=1/2.*[0 1 -1];
    otherplane(3,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[-1 1 1]/norm([-1 1 1]);
    splitdir(3,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    burg1(4,:)=1/6.*[1 2 1];
    burg2(4,:)=1/2.*[-1 0 -1];
    otherplane(4,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[-1 1 1]/norm([-1 1 1]);
    splitdir(4,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    burg1(5,:)=1/6.*[1 -1 -2];
    burg2(5,:)=1/2.*[-1 1 0];
    otherplane(5,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 1 -1]/norm([1 1 -1]);
    splitdir(5,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    if(norm(cross(tangent,([1 1 0]/(norm([1 1 0])))))<rntol)
        burg1(1,:)=1/6.*[-1 1 -2];
        burg2(1,:)=1/6.*[-1 1 0];
        otherplane(1,:)=[-1 1 1]/norm([-1 1 1]);
        splitdir(1,:)=[-1 1 -2]/norm([-1 1 -2]);
        burg1(2,:)=1/6.*[-1 1 -2];
        burg2(2,:)=1/6.*[-1 1 0];
        otherplane(2,:)=[-1 1 1]/norm([-1 1 1]);
        splitdir(2,:)=[1 -1 2]/norm([1 -1 2]);
        
    elseif(norm(cross(tangent,([0 1 1]/(norm([0 1 1])))))<rntol)
        burg1(1,:)=1/6.*[-2 1 -1];
        burg2(1,:)=1/6.*[0 1 -1];
        otherplane(1,:)=[1 1 -1]/norm([1 1 -1]);
        splitdir(1,:)=[-2 1 -1]/norm([-2 1 -1]);
        burg1(2,:)=1/6.*[-2 1 -1];
        burg2(2,:)=1/6.*[0 1 -1];
        otherplane(2,:)=[1 1 -1]/norm([1 1 -1]);
        splitdir(2,:)=[2 -1 1]/norm([2 -1 1]);
        
    elseif(norm(cross(tangent,([-1 0 1]/(norm([-1 0 1])))))<rntol)
        burg1(1,:)=1/6.*[-1 2 -1];
        burg2(1,:)=1/6.*[-1 0 -1];
        otherplane(1,:)=[1 1 1]/norm([1 1 1]); 
        splitdir(1,:)=[-1 2 -1]/norm([-1 2 -1]);
        burg1(2,:)=1/6.*[-1 2 -1];
        burg2(2,:)=1/6.*[-1 0 -1];
        otherplane(2,:)=[1 1 1]/norm([1 1 1]); 
        splitdir(2,:)=[1 -2 1]/norm([1 -2 1]);
        
    end
elseif(norm(burgv-(1/3.*[1 -1 1]))<rntol)  %Aa frank partial
    burg1(3,:)=1/6.*[2 1 -1];
    burg2(3,:)=1/2.*[0 -1 1];
    otherplane(3,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[-1 1 1]/norm([-1 1 1]);
    splitdir(3,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    burg1(4,:)=1/6.*[-1 -2 -1];
    burg2(4,:)=1/2.*[1 0 1];
    otherplane(4,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[-1 1 1]/norm([-1 1 1]);
    splitdir(4,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    burg1(5,:)=1/6.*[-1 1 2];
    burg2(5,:)=1/2.*[1 -1 0];
    otherplane(5,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 1 -1]/norm([1 1 -1]);
    splitdir(5,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    if(norm(cross(tangent,([1 1 0]/(norm([1 1 0])))))<rntol)
        burg1(1,:)=1/6.*[1 -1 2];
        burg2(1,:)=1/6.*[1 -1 0];
        otherplane(1,:)=[-1 1 1]/norm([-1 1 1]);
        splitdir(1,:)=[-1 1 -2]/norm([-1 1 -2]);
        burg1(2,:)=1/6.*[1 -1 2];
        burg2(2,:)=1/6.*[1 -1 0];
        otherplane(2,:)=[-1 1 1]/norm([-1 1 1]);
        splitdir(2,:)=[1 -1 2]/norm([1 -1 2]);
               
    elseif(norm(cross(tangent,([0 1 1]/(norm([0 1 1])))))<rntol)
        burg1(1,:)=1/6.*[2 -1 1];
        burg2(1,:)=1/6.*[0 -1 1];
        otherplane(1,:)=[1 1 -1]/norm([1 1 -1]);
        splitdir(1,:)=[-2 1 -1]/norm([-2 1 -1]);
        burg1(2,:)=1/6.*[2 -1 1];
        burg2(2,:)=1/6.*[0 -1 1];
        otherplane(2,:)=[1 1 -1]/norm([1 1 -1]);
        splitdir(2,:)=[2 -1 1]/norm([2 -1 1]);
        
    elseif(norm(cross(tangent,([-1 0 1]/(norm([-1 0 1])))))<rntol)
        burg1(1,:)=1/6.*[1 -2 1];
        burg2(1,:)=1/6.*[1 0 1];
        otherplane(1,:)=[1 1 1]/norm([1 1 1]); 
        splitdir(1,:)=[-1 2 -1]/norm([-1 2 -1]);
        burg1(2,:)=1/6.*[1 -2 1];
        burg2(2,:)=1/6.*[1 0 1];
        otherplane(2,:)=[1 1 1]/norm([1 1 1]); 
        splitdir(2,:)=[1 -2 1]/norm([1 -2 1]);
        
    end
elseif(norm(burgv-(1/3.*[1 -1 -1]))<rntol)  %bB frank partial
    burg1(3,:)=1/6.*[-1 -2 1];
    burg2(3,:)=1/2.*[1 0 -1];
    otherplane(3,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 1 1]/norm([1 1 1]);
    splitdir(3,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
        
    burg1(4,:)=1/6.*[-1 1 -2];
    burg2(4,:)=1/2.*[1 -1 0];
    otherplane(4,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 1 1]/norm([1 1 1]);
    splitdir(4,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    burg1(5,:)=1/6.*[2 1 1];
    burg2(5,:)=1/2.*[0 -1 -1];
    otherplane(5,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 1 -1]/norm([1 1 -1]);
    splitdir(5,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    if(norm(cross(tangent,([0 -1 1]/(norm([0 -1 1])))))<rntol)
        burg1(1,:)=1/6.*[2 -1 -1];
        burg2(1,:)=1/6.*[0 -1 -1];
        otherplane(1,:)=[1 1 1]/norm([1 1 1]);
        splitdir(1,:)=[2 -1 -1]/norm([2 -1 -1]);
        burg1(2,:)=1/6.*[2 -1 -1];
        burg2(2,:)=1/6.*[0 -1 -1];
        otherplane(2,:)=[1 1 1]/norm([1 1 1]);
        splitdir(2,:)=[-2 1 1]/norm([-2 1 1]);
        
    elseif(norm(cross(tangent,([1 0 1]/(norm([1 0 1])))))<rntol)
        burg1(1,:)=1/6.*[1 -2 -1];
        burg2(1,:)=1/6.*[1 0 -1];
        otherplane(1,:)=[1 1 -1]/norm([1 1 -1]);
        splitdir(1,:)=[1 -2 -1]/norm([1 -2 -1]);
        burg1(2,:)=1/6.*[1 -2 -1];
        burg2(2,:)=1/6.*[1 0 -1];
        otherplane(2,:)=[1 1 -1]/norm([1 1 -1]);
        splitdir(2,:)=[-1 2 1]/norm([-1 2 1]);
                
    elseif(norm(cross(tangent,([1 1 0]/(norm([1 1 0])))))<rntol)
        burg1(1,:)=1/6.*[1 -1 -2];
        burg2(1,:)=1/6.*[1 -1 0];
        otherplane(1,:)=[1 -1 1]/norm([1 -1 1]); 
        splitdir(1,:)=[1 -1 -2]/norm([1 -1 -2]);
        burg1(2,:)=1/6.*[1 -1 -2];
        burg2(2,:)=1/6.*[1 -1 0];
        otherplane(2,:)=[1 -1 1]/norm([1 -1 1]); 
        splitdir(2,:)=[-1 1 2]/norm([-1 1 2]);
        
    end
elseif(norm(burgv-(1/3.*[-1 1 1]))<rntol)  %Bb frank partial
    burg1(3,:)=1/6.*[1 2 -1];
    burg2(3,:)=1/2.*[-1 0 1];
    otherplane(3,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 1 1]/norm([1 1 1]);
    splitdir(3,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
        
    burg1(4,:)=1/6.*[1 -1 2];
    burg2(4,:)=1/2.*[-1 1 0];
    otherplane(4,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 1 1]/norm([1 1 1]);
    splitdir(4,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
        
    burg1(5,:)=1/6.*[-2 -1 -1];
    burg2(5,:)=1/2.*[0 1 1];
    otherplane(5,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 1 -1]/norm([1 1 -1]);
    splitdir(5,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    if(norm(cross(tangent,([0 -1 1]/(norm([0 -1 1])))))<rntol)
        burg1(1,:)=1/6.*[-2 1 1];
        burg2(1,:)=1/6.*[0 1 1];
        otherplane(1,:)=[1 1 1]/norm([1 1 1]);
        splitdir(1,:)=[2 -1 -1]/norm([2 -1 -1]);
        burg1(2,:)=1/6.*[-2 1 1];
        burg2(2,:)=1/6.*[0 1 1];
        otherplane(2,:)=[1 1 1]/norm([1 1 1]);
        splitdir(2,:)=[-2 1 1]/norm([-2 1 1]);
        
        
    elseif(norm(cross(tangent,([1 0 1]/(norm([1 0 1])))))<rntol)
        burg1(1,:)=1/6.*[-1 2 1];
        burg2(1,:)=1/6.*[-1 0 1];
        otherplane(1,:)=[1 1 -1]/norm([1 1 -1]);
        splitdir(1,:)=[1 -2 -1]/norm([1 -2 -1]);
        burg1(2,:)=1/6.*[-1 2 1];
        burg2(2,:)=1/6.*[-1 0 1];
        otherplane(2,:)=[1 1 -1]/norm([1 1 -1]);
        splitdir(2,:)=[-1 2 1]/norm([-1 2 1]);
        
        
    elseif(norm(cross(tangent,([1 1 0]/(norm([1 1 0])))))<rntol)
        burg1(1,:)=1/6.*[-1 1 2];
        burg2(1,:)=1/6.*[-1 1 0];
        otherplane(1,:)=[1 -1 1]/norm([1 -1 1]); 
        splitdir(1,:)=[1 -1 -2]/norm([1 -1 -2]);
        burg1(2,:)=1/6.*[-1 1 2];
        burg2(2,:)=1/6.*[-1 1 0];
        otherplane(2,:)=[1 -1 1]/norm([1 -1 1]); 
        splitdir(2,:)=[-1 1 2]/norm([-1 1 2]);
        
    end
elseif(norm(burgv-(1/3.*[1 1 -1]))<rntol)  %cC frank partial
    burg1(3,:)=1/6.*[2 -1 1];
    burg2(3,:)=1/2.*[0 1 -1];
    otherplane(3,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 1 1]/norm([1 1 1]);
    splitdir(3,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
        
    burg1(4,:)=1/6.*[-1 2 1];
    burg2(4,:)=1/2.*[1 0 -1];
    otherplane(4,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 1 1]/norm([1 1 1]);
    splitdir(4,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    burg1(5,:)=1/6.*[-1 -1 -2];
    burg2(5,:)=1/2.*[1 1 0];
    otherplane(5,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 -1 1]/norm([1 -1 1]);
    splitdir(5,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    if(norm(cross(tangent,([1 -1 0]/(norm([1 -1 0])))))<rntol)
        burg1(1,:)=1/6.*[1 1 -2];
        burg2(1,:)=1/6.*[1 1 0];
        otherplane(1,:)=[1 1 1]/norm([1 1 1]);
        splitdir(1,:)=[-1 -1 2]/norm([-1 -1 2]);
        burg1(2,:)=1/6.*[1 1 -2];
        burg2(2,:)=1/6.*[1 1 0];
        otherplane(2,:)=[1 1 1]/norm([1 1 1]);
        splitdir(2,:)=[1 1 -2]/norm([1 1 -2]);
        
        
    elseif(norm(cross(tangent,([0 1 1]/(norm([0 1 1])))))<rntol)
        burg1(1,:)=1/6.*[2 1 -1];
        burg2(1,:)=1/6.*[0 1 -1];
        otherplane(1,:)=[1 -1 1]/norm([1 -1 1]);
        splitdir(1,:)=[-2 -1 1]/norm([-2 -1 1]);
        burg1(2,:)=1/6.*[2 1 -1];
        burg2(2,:)=1/6.*[0 1 -1];
        otherplane(2,:)=[1 -1 1]/norm([1 -1 1]);
        splitdir(2,:)=[2 1 -1]/norm([2 1 -1]);
        
        
    elseif(norm(cross(tangent,([1 0 1]/(norm([1 0 1])))))<rntol)
        burg1(1,:)=1/6.*[1 2 -1];
        burg2(1,:)=1/6.*[1 0 -1];
        otherplane(1,:)=[-1 1 1]/norm([-1 1 1]);
        splitdir(1,:)=[-1 -2 1]/norm([-1 -2 1]);
        burg1(2,:)=1/6.*[1 2 -1];
        burg2(2,:)=1/6.*[1 0 -1];
        otherplane(2,:)=[-1 1 1]/norm([-1 1 1]);
        splitdir(2,:)=[1 2 -1]/norm([1 2 -1]);
        
        
    end
elseif(norm(burgv-(1/3.*[-1 -1 1]))<rntol)  %Cc frank partial
    burg1(3,:)=1/6.*[-2 1 -1];
    burg2(3,:)=1/2.*[0 -1 1];
    otherplane(3,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 1 1]/norm([1 1 1]);
    splitdir(3,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
        
    burg1(4,:)=1/6.*[1 -2 -1];
    burg2(4,:)=1/2.*[-1 0 1];
    otherplane(4,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 1 1]/norm([1 1 1]);
    splitdir(4,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    burg1(5,:)=1/6.*[1 1 2];
    burg2(5,:)=1/2.*[-1 -1 0];
    otherplane(5,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 -1 1]/norm([1 -1 1]);
    splitdir(5,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
        
    if(norm(cross(tangent,([1 -1 0]/(norm([1 -1 0])))))<rntol)
        burg1(1,:)=1/6.*[-1 -1 2];
        burg2(1,:)=1/6.*[-1 -1 0];
        otherplane(1,:)=[1 1 1]/norm([1 1 1]);
        splitdir(1,:)=[-1 -1 2]/norm([-1 -1 2]);
        burg1(2,:)=1/6.*[-1 -1 2];
        burg2(2,:)=1/6.*[-1 -1 0];
        otherplane(2,:)=[1 1 1]/norm([1 1 1]);
        splitdir(2,:)=[1 1 -2]/norm([1 1 -2]);
            
    elseif(norm(cross(tangent,([0 1 1]/(norm([0 1 1])))))<rntol)
        burg1(1,:)=1/6.*[-2 -1 1];
        burg2(1,:)=1/6.*[0 -1 1];
        otherplane(1,:)=[1 -1 1]/norm([1 -1 1]);
        splitdir(1,:)=[-2 -1 1]/norm([-2 -1 1]);
        burg1(2,:)=1/6.*[-2 -1 1];
        burg2(2,:)=1/6.*[0 -1 1];
        otherplane(2,:)=[1 -1 1]/norm([1 -1 1]);
        splitdir(2,:)=[2 1 -1]/norm([2 1 -1]);
        
        
    elseif(norm(cross(tangent,([1 0 1]/(norm([1 0 1])))))<rntol)
        burg1(1,:)=1/6.*[-1 -2 1];
        burg2(1,:)=1/6.*[-1 0 1];
        otherplane(1,:)=[-1 1 1]/norm([-1 1 1]);
        splitdir(1,:)=[-1 -2 1]/norm([-1 -2 1]);
        burg1(2,:)=1/6.*[-1 -2 1];
        burg2(2,:)=1/6.*[-1 0 1];
        otherplane(2,:)=[-1 1 1]/norm([-1 1 1]);
        splitdir(2,:)=[1 2 -1]/norm([1 2 -1]);
        
    end
elseif(norm(burgv-(1/3.*[1 1 1]))<rntol)  %dD frank partial
    burg1(3,:)=1/6.*[-1 2 -1];
    burg2(3,:)=1/2.*[1 0 1];
    otherplane(3,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[-1 1 1]/norm([-1 1 1]);
    splitdir(3,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
        
    burg1(4,:)=1/6.*[-1 -1 2];
    burg2(4,:)=1/2.*[1 1 0];
    otherplane(4,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[-1 1 1]/norm([-1 1 1]);
    splitdir(4,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    burg1(5,:)=1/6.*[2 -1 -1];
    burg2(5,:)=1/2.*[0 1 1];
    otherplane(5,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 -1 1]/norm([1 -1 1]);
    splitdir(5,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
    if(norm(cross(tangent,([0 -1 1]/(norm([0 -1 1])))))<rntol)
        burg1(1,:)=1/6.*[2 1 1];
        burg2(1,:)=1/6.*[0 1 1];
        otherplane(1,:)=[-1 1 1]/norm([-1 1 1]);
        splitdir(1,:)=[2 1 1]/norm([2 1 1]);
        burg1(2,:)=1/6.*[2 1 1];
        burg2(2,:)=1/6.*[0 1 1];
        otherplane(2,:)=[-1 1 1]/norm([-1 1 1]);
        splitdir(2,:)=[-2 -1 -1]/norm([-2 -1 -1]);
        
        
    elseif(norm(cross(tangent,([-1 0 1]/(norm([-1 0 1])))))<rntol)
        burg1(1,:)=1/6.*[1 2 1];
        burg2(1,:)=1/6.*[1 0 1];
        otherplane(1,:)=[1 -1 1]/norm([1 -1 1]);
        splitdir(1,:)=[1 2 1]/norm([1 2 1]);
        burg1(2,:)=1/6.*[1 2 1];
        burg2(2,:)=1/6.*[1 0 1];
        otherplane(2,:)=[1 -1 1]/norm([1 -1 1]);
        splitdir(2,:)=[-1 -2 -1]/norm([-1 -2 -1]);
        
        
    elseif(norm(cross(tangent,([1 -1 0]/(norm([1 -1 0])))))<rntol)
        burg1(1,:)=1/6.*[1 1 2];
        burg2(1,:)=1/6.*[1 1 0];
        otherplane(1,:)=[1 1 -1]/norm([1 1 -1]); 
        splitdir(1,:)=[1 1 2]/norm([1 1 2]);
        burg1(2,:)=1/6.*[1 1 2];
        burg2(2,:)=1/6.*[1 1 0];
        otherplane(2,:)=[1 1 -1]/norm([1 1 -1]); 
        splitdir(2,:)=[-1 -1 -2]/norm([-1 -1 -2]);
        
        
    end
   
elseif(norm(burgv-(1/3.*[-1 -1 -1]))<rntol)  %Dd frank partial
    burg1(3,:)=1/6.*[1 -2 1];
    burg2(3,:)=1/2.*[-1 0 -1];
    otherplane(3,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[-1 1 1]/norm([-1 1 1]);
    splitdir(3,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
        
    burg1(4,:)=1/6.*[1 1 -2];
    burg2(4,:)=1/2.*[-1 -1 0];
    otherplane(4,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[-1 1 1]/norm([-1 1 1]);
    splitdir(4,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
        
    burg1(5,:)=1/6.*[-2 1 1];
    burg2(5,:)=1/2.*[0 -1 -1];
    otherplane(5,:)=cross(tangent,burg2(3,:))/norm(cross(tangent,burg2(3,:)));%[1 -1 1]/norm([1 -1 1]);
    splitdir(5,:)=cross(tangent,stackforcevec)/norm(cross(tangent,stackforcevec));
        
    if(norm(cross(tangent,([0 -1 1]/(norm([0 -1 1])))))<rntol)
        burg1(1,:)=1/6.*[-2 -1 -1];
        burg2(1,:)=1/6.*[0 -1 -1];
        otherplane(1,:)=[-1 1 1]/norm([-1 1 1]);
        splitdir(1,:)=[2 1 1]/norm([2 1 1]);
        burg1(2,:)=1/6.*[-2 -1 -1];
        burg2(2,:)=1/6.*[0 -1 -1];
        otherplane(2,:)=[-1 1 1]/norm([-1 1 1]);
        splitdir(2,:)=[-2 -1 -1]/norm([-2 -1 -1]);
        
        
    elseif(norm(cross(tangent,([-1 0 1]/(norm([-1 0 1])))))<rntol)
        burg1(1,:)=1/6.*[-1 -2 -1];
        burg2(1,:)=1/6.*[-1 0 -1];
        otherplane(1,:)=[1 -1 1]/norm([1 -1 1]);
        splitdir(1,:)=[1 2 1]/norm([1 2 1]);
        burg1(2,:)=1/6.*[-1 -2 -1];
        burg2(2,:)=1/6.*[-1 0 -1];
        otherplane(2,:)=[1 -1 1]/norm([1 -1 1]);
        splitdir(2,:)=[-1 -2 -1]/norm([-1 -2 -1]);
        
        
    elseif(norm(cross(tangent,([1 -1 0]/(norm([1 -1 0])))))<rntol)
        burg1(1,:)=1/6.*[-1 -1 -2];
        burg2(1,:)=1/6.*[-1 -1 0];
        otherplane(1,:)=[1 1 -1]/norm([1 1 -1]); 
        splitdir(1,:)=[1 1 2]/norm([1 1 2]);
        burg1(2,:)=1/6.*[-1 -1 -2];
        burg2(2,:)=1/6.*[-1 -1 0];
        otherplane(2,:)=[1 1 -1]/norm([1 1 -1]); 
        splitdir(2,:)=[-1 -1 -2]/norm([-1 -1 -2]);
        
        
    end
end


    
    
        
    