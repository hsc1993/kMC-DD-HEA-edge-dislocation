%This function returns the partial Burgers vectors depending on the perfect Burgers vector and the glide plane

function [burgpart1,burgpart2,burgpartcsp1,burgpartcsp2,splitdirgp,splitdircsp,crossslipplane]=burgervecpartials(burgiperf,tangperfbt,normplaneperf)
tangperf=tangperfbt(1,:)/norm(tangperfbt(1,:));
burgperf=burgiperf(1,:);
burgpart1=zeros(1,3);
burgpart2=zeros(1,3);
burgpartcsp1=zeros(1,3);
burgpartcsp2=zeros(1,3);
splitdirgp=zeros(1,3);
splitdircsp=zeros(1,3);
crossslipplane=zeros(1,3);
eps=1e-1;

normal1=cross(burgiperf(1,:),tangperfbt(1,:));
normal2=cross(burgiperf(2,:),tangperfbt(2,:));
if(norm(normal1)>eps)
    normal=normal1/(norm(normal1));
elseif(norm(normal2)>eps)
    normal=normal2/(norm(normal2));
else
    normal=normplaneperf(1,:)/norm(normplaneperf(1,:));
end
%Now we have to introduce all the possibilities
%We fix the normal and change the Burgers vectors
%This is for the [-1 -1 -1] glide plane OR [1 1 1]
if((norm(cross(normal,([-1 -1 -1]/(norm([-1 -1 -1])))))<eps)|(norm(cross(normal,([1 1 1]/(norm([1 1 1])))))<eps))
    if(burgperf==(1/2.*[1 -1 0]))
        burgpart1=1/6.*[1 -2 1];
        burgpart2=1/6.*[2 -1 -1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[2 -1 1];
            burgpartcsp2=1/6.*[1 -2 -1];
            crossslipplane=[1 1 -1]/norm([1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 -1 1]/norm([0 -1 1])))<eps)
            burgpartcsp1=1/6.*[1 -1 2];
            burgpartcsp2=1/3.*[1 -1 -1];
            crossslipplane=[-1 1 1]/norm([-1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[-1 0 1]/norm([-1 0 1])))<eps)
            burgpartcsp1=1/6.*[1 -1 -2];
            burgpartcsp2=1/3.*[1 -1 1];
            crossslipplane=[1 -1 1]/norm([1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end
        
    elseif(burgperf==(1/2.*[-1 1 0]))
        burgpart1=1/6.*[-2 1 1];
        burgpart2=1/6.*[-1 2 -1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[-2 1 -1];
            burgpartcsp2=1/6.*[-1 2 1];
            crossslipplane=[1 1 -1]/norm([1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 -1 1]/norm([0 -1 1])))<eps)
            burgpartcsp1=1/6.*[-1 1 -2];
            burgpartcsp2=1/3.*[-1 1 1];
            crossslipplane=[1 -1 -1]/norm([1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 0 -1]/norm([1 0 -1])))<eps)
            burgpartcsp1=1/6.*[-1 1 2];
            burgpartcsp2=1/3.*[-1 1 -1];
            crossslipplane=[-1 1 -1]/norm([-1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end
        
    elseif(burgperf==(1/2.*[-1 0 1]))
        burgpart1=1/6.*[-2 1 1];
        burgpart2=1/6.*[-1 -1 2];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[-2 -1 1];
            burgpartcsp2=1/6.*[-1 1 2];
            crossslipplane=[1 -1 1]/norm([1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 -1 1]/norm([0 -1 1])))<eps)
            burgpartcsp1=1/6.*[-1 -2 1];
            burgpartcsp2=1/3.*[-1 1 1];
            crossslipplane=[1 -1 -1]/norm([1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 -1 0]/norm([1 -1 0])))<eps)
            burgpartcsp1=1/6.*[-1 2 1];
            burgpartcsp2=1/3.*[-1 -1 1];
            crossslipplane=[-1 -1 1]/norm([-1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end
        
        
    elseif(burgperf==(1/2.*[1 0 -1]))
        burgpart1=1/6.*[1 1 -2];
        burgpart2=1/6.*[2 -1 -1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[2 1 -1];
            burgpartcsp2=1/6.*[1 -1 -2];
            crossslipplane=[1 -1 1]/norm([1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 -1 1]/norm([0 -1 1])))<eps)
            burgpartcsp1=1/6.*[1 2 -1];
            burgpartcsp2=1/3.*[1 -1 -1];
            crossslipplane=[1 -1 -1]/norm([1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 -1 0]/norm([1 -1 0])))<eps)
            burgpartcsp1=1/6.*[1 -2 -1];
            burgpartcsp2=1/3.*[1 1 -1];
            crossslipplane=[-1 -1 1]/norm([-1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
    elseif(burgperf==(1/2.*[0 -1 1]))
        burgpart1=1/6.*[1 -2 1];
        burgpart2=1/6.*[-1 -1 2];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[1 -1 2];
            burgpartcsp2=1/6.*[-1 -2 1];
            crossslipplane=[-1 1 1]/norm([-1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[-1 0 1]/norm([-1 0 1])))<eps)
            burgpartcsp1=1/6.*[-2 -1 1];
            burgpartcsp2=1/3.*[1 -1 1];
            crossslipplane=[1 -1 1]/norm([1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 -1 0]/norm([1 -1 0])))<eps)
            burgpartcsp1=1/6.*[2 -1 1];
            burgpartcsp2=1/3.*[-1 -1 1];
            crossslipplane=[-1 -1 1]/norm([-1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
    elseif(burgperf==(1/2.*[0 1 -1]))
        burgpart1=1/6.*[1 1 -2];
        burgpart2=1/6.*[-1 2 -1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[-1 1 -2];
            burgpartcsp2=1/6.*[1 2 -1];
            crossslipplane=[-1 1 1]/norm([-1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[-1 0 1]/norm([-1 0 1])))<eps)
            burgpartcsp1=1/6.*[2 1 -1];
            burgpartcsp2=1/3.*[-1 1 -1];
            crossslipplane=[-1 1 -1]/norm([-1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 -1 0]/norm([1 -1 0])))<eps)
            burgpartcsp1=1/6.*[-2 1 -1];
            burgpartcsp2=1/3.*[1 1 -1];
            crossslipplane=[1 1 -1]/norm([1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
    end
    
    %This is for the plane [1 -1 1] OR [-1 1 -1]
elseif((norm(cross(normal,([1 -1 1]/(norm([1 -1 1])))))<eps)|(norm(cross(normal,([-1 1 -1]/(norm([-1 1 -1])))))<eps))
    
    if(burgperf==(1/2.*[1 1 0]))
        burgpart1=1/6.*[2 1 -1];
        burgpart2=1/6.*[1 2 1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[1 2 -1];
            burgpartcsp2=1/6.*[2 1 1];
            crossslipplane=[-1 1 1]/norm([-1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[-1 0 1]/norm([-1 0 1])))<eps)
            burgpartcsp1=1/6.*[1 1 -2];
            burgpartcsp2=1/3.*[1 1 1];
            crossslipplane=[1 1 1]/norm([1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 1 1]/norm([0 1 1])))<eps)
            burgpartcsp1=1/6.*[1 1 2];
            burgpartcsp2=1/3.*[-1 -1 1];
            crossslipplane=[-1 -1 1]/norm([-1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end
        
    elseif(burgperf==(1/2.*[-1 -1 0]))
        burgpart1=1/6.*[-1 -2 -1];
        burgpart2=1/6.*[-2 -1 1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[-1 -2 1];
            burgpartcsp2=1/6.*[-2 -1 -1];
            crossslipplane=[-1 1 1]/norm([-1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[-1 0 1]/norm([-1 0 1])))<eps)
            burgpartcsp1=1/6.*[-1 -1 2];
            burgpartcsp2=1/3.*[-1 -1 -1];
            crossslipplane=[-1 -1 -1]/norm([-1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 1 1]/norm([0 1 1])))<eps)
            burgpartcsp1=1/6.*[-1 -1 -2];
            burgpartcsp2=1/3.*[1 1 -1];
            crossslipplane=[1 1 -1]/norm([1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
        
    elseif(burgperf==(1/2.*[-1 0 1]))
        burgpart1=1/6.*[-1 1 2];
        burgpart2=1/6.*[-2 -1 1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[-2 1 1];
            burgpartcsp2=1/6.*[-1 -1 2];
            crossslipplane=[-1 -1 -1]/norm([-1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 1 0]/norm([1 1 0])))<eps)
            burgpartcsp1=1/6.*[-1 -2 1];
            burgpartcsp2=1/3.*[-1 1 1];
            crossslipplane=[-1 1 1]/norm([-1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 1 1]/norm([0 1 1])))<eps)
            burgpartcsp1=1/6.*[-1 2 1];
            burgpartcsp2=1/3.*[-1 -1 1];
            crossslipplane=[-1 -1 1]/norm([-1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
        
    elseif(burgperf==(1/2.*[1 0 -1]))
        burgpart1=1/6.*[2 1 -1];
        burgpart2=1/6.*[1 -1 -2];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[2 -1 -1];
            burgpartcsp2=1/6.*[1 1 -2];
            crossslipplane=[-1 -1 -1]/norm([-1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 1 0]/norm([1 1 0])))<eps)
            burgpartcsp1=1/6.*[1 2 -1];
            burgpartcsp2=1/3.*[1 -1 -1];
            crossslipplane=[1 -1 -1]/norm([1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 1 1]/norm([0 1 1])))<eps)
            burgpartcsp1=1/6.*[1 -2 -1];
            burgpartcsp2=1/3.*[1 1 -1];
            crossslipplane=[1 1 -1]/norm([1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
        
    elseif(burgperf==(1/2.*[0 1 1]))
        burgpart1=1/6.*[-1 1 2];
        burgpart2=1/6.*[1 2 1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[-1 2 1];
            burgpartcsp2=1/6.*[1 1 2];
            crossslipplane=[1 1 -1]/norm([1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 1 0]/norm([1 1 0])))<eps)
            burgpartcsp1=1/6.*[2 1 1];
            burgpartcsp2=1/3.*[-1 1 1];
            crossslipplane=[-1 1 1]/norm([-1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[-1 0 1]/norm([-1 0 1])))<eps)
            burgpartcsp1=1/6.*[-2 1 1];
            burgpartcsp2=1/3.*[1 1 1];
            crossslipplane=[1 1 1]/norm([1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
        
    elseif(burgperf==(1/2.*[0 -1 -1]))
        burgpart1=1/6.*[-1 -2 -1];
        burgpart2=1/6.*[1 -1 -2];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[1 -2 -1];
            burgpartcsp2=1/6.*[-1 -1 -2];
            crossslipplane=[1 1 -1]/norm([1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 1 0]/norm([1 1 0])))<eps)
            burgpartcsp1=1/6.*[-2 -1 -1];
            burgpartcsp2=1/3.*[1 -1 -1];
            crossslipplane=[1 -1 -1]/norm([1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[-1 0 1]/norm([-1 0 1])))<eps)
            burgpartcsp1=1/6.*[2 -1 -1];
            burgpartcsp2=1/3.*[-1 -1 -1];
            crossslipplane=[-1 -1 -1]/norm([-1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
    end
    
    %This is for the plane [1 1 -1] OR [-1 -1 1]
elseif(norm(cross(normal,([1 1 -1]/(norm([1 1 -1])))))<eps|(norm(cross(normal,([-1 -1 1]/(norm([-1 -1 1])))))<eps))
    
    if(burgperf==(1/2.*[1 0 1]))
        burgpart1=1/6.*[2 -1 1];
        burgpart2=1/6.*[1 1 2];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[1 -1 2];
            burgpartcsp2=1/6.*[2 1 1];
            crossslipplane=[-1 1 1]/norm([-1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 -1 0]/norm([1 -1 0])))<eps)
            burgpartcsp1=1/6.*[1 -2 1];
            burgpartcsp2=1/3.*[1 1 1];
            crossslipplane=[1 1 1]/norm([1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 1 1]/norm([0 1 1])))<eps)
            burgpartcsp1=1/6.*[1 2 1];
            burgpartcsp2=1/3.*[1 -1 1];
            crossslipplane=[1 -1 1]/norm([1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
        
    elseif(burgperf==(1/2.*[-1 0 -1]))
        burgpart1=1/6.*[-1 -1 -2];
        burgpart2=1/6.*[-2 1 -1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[-1 1 -2];
            burgpartcsp2=1/6.*[-2 -1 -1];
            crossslipplane=[-1 1 1]/norm([-1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 -1 0]/norm([1 -1 0])))<eps)
            burgpartcsp1=1/6.*[-1 2 -1];
            burgpartcsp2=1/3.*[-1 -1 -1];
            crossslipplane=[-1 -1 -1]/norm([-1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 1 1]/norm([0 1 1])))<eps)
            burgpartcsp1=1/6.*[-1 -2 -1];
            burgpartcsp2=1/3.*[-1 1 -1];
            crossslipplane=[-1 1 -1]/norm([-1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
        
    elseif(burgperf==(1/2.*[0 1 1]))
        burgpart1=1/6.*[-1 2 1];
        burgpart2=1/6.*[1 1 2];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[-1 1 2];
            burgpartcsp2=1/6.*[1 2 1];
            crossslipplane=[1 -1 1]/norm([1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 -1 0]/norm([1 -1 0])))<eps)
            burgpartcsp1=1/6.*[-2 1 1];
            burgpartcsp2=1/3.*[1 1 1];
            crossslipplane=[1 1 1]/norm([1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 0 1]/norm([1 0 1])))<eps)
            burgpartcsp1=1/6.*[2 1 1];
            burgpartcsp2=1/3.*[-1 1 1];
            crossslipplane=[-1 1 1]/norm([-1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
        
    elseif(burgperf==(1/2.*[0 -1 -1]))
        burgpart1=1/6.*[-1 -1 -2];
        burgpart2=1/6.*[1 -2 -1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[1 -1 -2];
            burgpartcsp2=1/6.*[-1 -2 -1];
            crossslipplane=[1 -1 1]/norm([1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 -1 0]/norm([1 -1 0])))<eps)
            burgpartcsp1=1/6.*[2 -1 -1];
            burgpartcsp2=1/3.*[-1 -1 -1];
            crossslipplane=[-1 -1 -1]/norm([-1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 0 1]/norm([1 0 1])))<eps)
            burgpartcsp1=1/6.*[-2 -1 -1];
            burgpartcsp2=1/3.*[1 -1 -1];
            crossslipplane=[1 -1 -1]/norm([1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
        
    elseif(burgperf==(1/2.*[1 -1 0]))
        burgpart1=1/6.*[2 -1 1];
        burgpart2=1/6.*[1 -2 -1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[1 -2 1];
            burgpartcsp2=1/6.*[2 -1 -1];
            crossslipplane=[-1 -1 -1]/norm([-1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 1 1]/norm([0 1 1])))<eps)
            burgpartcsp1=1/6.*[1 -1 -2];
            burgpartcsp2=1/3.*[1 -1 1];
            crossslipplane=[1 -1 1]/norm([1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 0 1]/norm([1 0 1])))<eps)
            burgpartcsp1=1/6.*[1 -1 2];
            burgpartcsp2=1/3.*[1 -1 -1];
            crossslipplane=[1 -1 -1]/norm([1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
        
    elseif(burgperf==(1/2.*[-1 1 0]))
        burgpart1=1/6.*[-1 2 1];
        burgpart2=1/6.*[-2 1 -1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[-1 2 -1];
            burgpartcsp2=1/6.*[-2 1 1];
            crossslipplane=[-1 -1 -1]/norm([-1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 1 1]/norm([0 1 1])))<eps)
            burgpartcsp1=1/6.*[-1 1 2];
            burgpartcsp2=1/3.*[-1 1 -1];
            crossslipplane=[-1 1 -1]/norm([-1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 0 1]/norm([1 0 1])))<eps)
            burgpartcsp1=1/6.*[-1 1 -2];
            burgpartcsp2=1/3.*[-1 1 1];
            crossslipplane=[-1 1 1]/norm([-1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
    end
    
    %This is for the plane [-1 1 1] OR [1 -1 -1]
elseif((norm(cross(normal,([-1 1 1]/(norm([-1 1 1])))))<eps)|(norm(cross(normal,([1 -1 -1]/(norm([1 -1 -1])))))<eps))
    
    if(burgperf==(1/2.*[1 1 0]))
        burgpart1=1/6.*[1 2 -1];
        burgpart2=1/6.*[2 1 1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[2 1 -1];
            burgpartcsp2=1/6.*[1 2 1];
            crossslipplane=[1 -1 1]/norm([1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 -1 1]/norm([0 -1 1])))<eps)
            burgpartcsp1=1/6.*[1 1 -2];
            burgpartcsp2=1/3.*[1 1 1];
            crossslipplane=[1 1 1]/norm([1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 0 1]/norm([1 0 1])))<eps)
            burgpartcsp1=1/6.*[1 1 2];
            burgpartcsp2=1/3.*[1 1 -1];
            crossslipplane=[1 1 -1]/norm([1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
        
    elseif(burgperf==(1/2.*[-1 -1 0]))
        burgpart1=1/6.*[-2 -1 -1];
        burgpart2=1/6.*[-1 -2 1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[-2 -1 1];
            burgpartcsp2=1/6.*[-1 -2 -1];
            crossslipplane=[1 -1 1]/norm([1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 -1 1]/norm([0 -1 1])))<eps)
            burgpartcsp1=1/6.*[-1 -1 2];
            burgpartcsp2=1/3.*[-1 -1 -1];
            crossslipplane=[-1 -1 -1]/norm([-1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 0 1]/norm([1 0 1])))<eps)
            burgpartcsp1=1/6.*[-1 -1 -2];
            burgpartcsp2=1/3.*[-1 -1 1];
            crossslipplane=[-1 -1 1]/norm([-1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
        
    elseif(burgperf==(1/2.*[1 0 1]))
        burgpart1=1/6.*[1 -1 2];
        burgpart2=1/6.*[2 1 1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[2 -1 1];
            burgpartcsp2=1/6.*[1 1 2];
            crossslipplane=[1 1 -1]/norm([1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 -1 1]/norm([0 -1 1])))<eps)
            burgpartcsp1=1/6.*[1 -2 1];
            burgpartcsp2=1/3.*[1 1 1];
            crossslipplane=[1 1 1]/norm([1 1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 1 0]/norm([1 1 0])))<eps)
            burgpartcsp1=1/6.*[1 2 1];
            burgpartcsp2=1/3.*[1 -1 1];
            crossslipplane=[1 -1 1]/norm([1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
        
    elseif(burgperf==(1/2.*[-1 0 -1]))
        burgpart1=1/6.*[-2 -1 -1];
        burgpart2=1/6.*[-1 1 -2];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[-2 1 -1];
            burgpartcsp2=1/6.*[-1 -1 -2];
            crossslipplane=[1 1 -1]/norm([1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[0 -1 1]/norm([0 -1 1])))<eps)
            burgpartcsp1=1/6.*[-1 2 -1];
            burgpartcsp2=1/3.*[-1 -1 -1];
            crossslipplane=[-1 -1 -1]/norm([-1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 1 0]/norm([1 1 0])))<eps)
            burgpartcsp1=1/6.*[-1 -2 -1];
            burgpartcsp2=1/3.*[-1 1 -1];
            crossslipplane=[-1 1 -1]/norm([-1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end
        
        
    elseif(burgperf==(1/2.*[0 -1 1]))
        burgpart1=1/6.*[1 -1 2];
        burgpart2=1/6.*[-1 -2 1];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[1 -2 1];
            burgpartcsp2=1/6.*[-1 -1 2];
            crossslipplane=[-1 -1 -1]/norm([-1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 1 0]/norm([1 1 0])))<eps)
            burgpartcsp1=1/6.*[-2 -1 1];
            burgpartcsp2=1/3.*[1 -1 1];
            crossslipplane=[1 -1 1]/norm([1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 0 1]/norm([1 0 1])))<eps)
            burgpartcsp1=1/6.*[2 -1 1];
            burgpartcsp2=1/3.*[-1 -1 1];
            crossslipplane=[-1 -1 1]/norm([-1 -1 1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
        
    elseif(burgperf==(1/2.*[0 1 -1]))
        burgpart1=1/6.*[1 2 -1];
        burgpart2=1/6.*[-1 1 -2];
        splitdirgp=cross(tangperf,normal)/norm(cross(tangperf,normal));
        if(norm(cross(tangperf,burgperf))<eps)
            burgpartcsp1=1/6.*[-1 2 -1];
            burgpartcsp2=1/6.*[1 1 -2];
            crossslipplane=[-1 -1 -1]/norm([-1 -1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 1 0]/norm([1 1 0])))<eps)
            burgpartcsp1=1/6.*[2 1 -1];
            burgpartcsp2=1/3.*[-1 1 -1];
            crossslipplane=[-1 1 -1]/norm([-1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        elseif(norm(cross(tangperf,[1 0 1]/norm([1 0 1])))<eps)
            burgpartcsp1=1/6.*[-2 1 -1];
            burgpartcsp2=1/3.*[1 1 -1];
            crossslipplane=[1 1 -1]/norm([1 1 -1]);
            splitdircsp=cross(tangperf,crossslipplane)/norm(cross(tangperf,crossslipplane));
        end    
        
        
    end
end


