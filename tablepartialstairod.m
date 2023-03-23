%This function assigns the new Burgers vectors and add a new node

function [newburpar,newstairod,splitdir]=tablepartialstairod(burgv,tangent,rntol)

newburpar=zeros(1,3);
newstairod=zeros(1,3);
splitdir=zeros(1,3);

if((norm(cross(tangent/norm(tangent),[-1 0 1]/norm([-1 0 1])))<rntol)|(norm(cross(tangent/norm(tangent),[1 0 -1]/norm([1 0 -1])))<rntol))
    if(norm(burgv-1/6.*[-2 1 1])<rntol)
        newburpar=1/6.*[-1 1 2];
        newstairod=1/6.*[-1 0 -1];
        splitdir=cross([1 -1 1]/norm([1 -1 1]),tangent)/norm(cross([1 -1 1]/norm([1 -1 1]),tangent));
        
            
        
    elseif(norm(burgv-1/6.*[2 -1 -1])<rntol)
        newburpar=1/6.*[1 -1 -2];
        newstairod=1/6.*[1 0 1];
        splitdir=cross([1 -1 1]/norm([1 -1 1]),tangent)/norm(cross([1 -1 1]/norm([1 -1 1]),tangent));
        
    elseif(norm(burgv-1/6.*[-1 1 2])<rntol)
        newburpar=1/6.*[-2 1 1];
        newstairod=1/6.*[1 0 1];
        splitdir=cross([1 1 1]/norm([1 1 1]),tangent)/norm(cross([1 1 1]/norm([1 1 1]),tangent));
        
    elseif(norm(burgv-1/6.*[1 -1 -2])<rntol)
        newburpar=1/6.*[2 -1 -1];
        newstairod=1/6.*[-1 0 -1];
        splitdir=cross([1 1 1]/norm([1 1 1]),tangent)/norm(cross([1 1 1]/norm([1 1 1]),tangent));
        
    elseif(norm(burgv-1/6.*[1 1 -2])<rntol)
        newburpar=1/6.*[2 1 -1];
        newstairod=1/6.*[-1 0 -1];
        splitdir=cross([1 -1 1]/norm([1 -1 1]),tangent)/norm(cross([1 -1 1]/norm([1 -1 1]),tangent));
        
    elseif(norm(burgv-1/6.*[-1 -1 2])<rntol)
        newburpar=1/6.*[-2 -1 1];
        newstairod=1/6.*[1 0 1];
        splitdir=cross([1 -1 1]/norm([1 -1 1]),tangent)/norm(cross([1 -1 1]/norm([1 -1 1]),tangent));
        
    elseif(norm(burgv-1/6.*[2 1 -1])<rntol)
        newburpar=1/6.*[1 1 -2];
        newstairod=1/6.*[1 0 1];
        splitdir=cross([1 1 1]/norm([1 1 1]),tangent)/norm(cross([1 1 1]/norm([1 1 1]),tangent));  
        
    elseif(norm(burgv-1/6.*[-2 -1 1])<rntol)
        newburpar=1/6.*[-1 -1 2];
        newstairod=1/6.*[-1 0 -1];
        splitdir=cross([1 1 1]/norm([1 1 1]),tangent)/norm(cross([1 1 1]/norm([1 1 1]),tangent)); 
        
    end
    
elseif((norm(cross(tangent/norm(tangent),[1 -1 0]/norm([1 -1 0])))<rntol)|(norm(cross(tangent/norm(tangent),[-1 1 0]/norm([-1 1 0])))<rntol))
    if(norm(burgv-1/6.*[-2 1 1])<rntol)
        newburpar=1/6.*[-1 2 1];
        newstairod=1/6.*[-1 -1 0];
        splitdir=cross([1 1 -1]/norm([1 1 -1]),tangent)/norm(cross([1 1 -1]/norm([1 1 -1]),tangent));
        
    elseif(norm(burgv-1/6.*[2 -1 -1])<rntol)
        newburpar=1/6.*[1 -2 -1];
        newstairod=1/6.*[1 1 0];
        splitdir=cross([1 1 -1]/norm([1 1 -1]),tangent)/norm(cross([1 1 -1]/norm([1 1 -1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[-1 2 1])<rntol)
        newburpar=1/6.*[-2 1 1];
        newstairod=1/6.*[1 1 0];
        splitdir=cross([1 1 1]/norm([1 1 1]),tangent)/norm(cross([1 1 1]/norm([1 1 1]),tangent));   
        
    elseif(norm(burgv-1/6.*[1 -2 -1])<rntol)
        newburpar=1/6.*[2 -1 -1];
        newstairod=1/6.*[-1 -1 0];
        splitdir=cross([1 1 1]/norm([1 1 1]),tangent)/norm(cross([1 1 1]/norm([1 1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[1 -2 1])<rntol)
        newburpar=1/6.*[2 -1 1];
        newstairod=1/6.*[-1 -1 0];
        splitdir=cross([1 1 -1]/norm([1 1 -1]),tangent)/norm(cross([1 1 -1]/norm([1 1 -1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[-1 2 -1])<rntol)
        newburpar=1/6.*[-2 1 -1];
        newstairod=1/6.*[1 1 0];
        splitdir=cross([1 1 -1]/norm([1 1 -1]),tangent)/norm(cross([1 1 -1]/norm([1 1 -1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[2 -1 1])<rntol)
        newburpar=1/6.*[1 -2 1];
        newstairod=1/6.*[1 1 0];
        splitdir=cross([1 1 1]/norm([1 1 1]),tangent)/norm(cross([1 1 1]/norm([1 1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[-2 1 -1])<rntol)
        newburpar=1/6.*[-1 2 -1];
        newstairod=1/6.*[-1 -1 0];
        splitdir=cross([1 1 1]/norm([1 1 1]),tangent)/norm(cross([1 1 1]/norm([1 1 1]),tangent));  
        
    end
 elseif((norm(cross(tangent/norm(tangent),[0 -1 1]/norm([0 -1 1])))<rntol)|(norm(cross(tangent/norm(tangent),[0 1 -1]/norm([0 1 -1])))<rntol))
    if(norm(burgv-1/6.*[1 -2 1])<rntol)
        newburpar=1/6.*[1 -1 2];
        newstairod=1/6.*[0 -1 -1];
        splitdir=cross([-1 1 1]/norm([-1 1 1]),tangent)/norm(cross([-1 1 1]/norm([-1 1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[-1 2 -1])<rntol)
        newburpar=1/6.*[-1 1 -2];
        newstairod=1/6.*[0 1 1];
        splitdir=cross([-1 1 1]/norm([-1 1 1]),tangent)/norm(cross([-1 1 1]/norm([-1 1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[1 -1 2])<rntol)
        newburpar=1/6.*[1 -2 1];
        newstairod=1/6.*[0 1 1];
        splitdir=cross([1 1 1]/norm([1 1 1]),tangent)/norm(cross([1 1 1]/norm([1 1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[-1 1 -2])<rntol)
        newburpar=1/6.*[-1 2 -1];
        newstairod=1/6.*[0 -1 -1];
        splitdir=cross([1 1 1]/norm([1 1 1]),tangent)/norm(cross([1 1 1]/norm([1 1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[1 1 -2])<rntol)
        newburpar=1/6.*[1 2 -1];
        newstairod=1/6.*[0 -1 -1];
        splitdir=cross([-1 1 1]/norm([-1 1 1]),tangent)/norm(cross([-1 1 1]/norm([-1 1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[-1 -1 2])<rntol)
        newburpar=1/6.*[-1 -2 1];
        newstairod=1/6.*[0 1 1];
        splitdir=cross([-1 1 1]/norm([-1 1 1]),tangent)/norm(cross([-1 1 1]/norm([-1 1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[1 2 -1])<rntol)
        newburpar=1/6.*[1 1 -2];
        newstairod=1/6.*[0 1 1];
        splitdir=cross([1 1 1]/norm([1 1 1]),tangent)/norm(cross([1 1 1]/norm([1 1 1]),tangent));
        
    elseif(norm(burgv-1/6.*[-1 -2 1])<rntol)
        newburpar=1/6.*[-1 -1 2];
        newstairod=1/6.*[0 -1 -1];
        splitdir=cross([1 1 1]/norm([1 1 1]),tangent)/norm(cross([1 1 1]/norm([1 1 1]),tangent)); 
        
    end
           
 elseif((norm(cross(tangent/norm(tangent),[1 1 0]/norm([1 1 0])))<rntol)|(norm(cross(tangent/norm(tangent),[-1 -1 0]/norm([-1 -1 0])))<rntol))
    if(norm(burgv-1/6.*[1 2 -1])<rntol)
        newburpar=1/6.*[2 1 -1];
        newstairod=1/6.*[-1 1 0];
        splitdir=cross([1 -1 1]/norm([1 -1 1]),tangent)/norm(cross([1 -1 1]/norm([1 -1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[-1 -2 1])<rntol)
        newburpar=1/6.*[-2 -1 1];
        newstairod=1/6.*[1 -1 0];
        splitdir=cross([1 -1 1]/norm([1 -1 1]),tangent)/norm(cross([1 -1 1]/norm([1 -1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[2 1 -1])<rntol)
        newburpar=1/6.*[1 2 -1];
        newstairod=1/6.*[1 -1 0];
        splitdir=cross([-1 1 1]/norm([-1 1 1]),tangent)/norm(cross([-1 1 1]/norm([-1 1 1]),tangent));  
        
    elseif(norm(burgv-1/6.*[-2 -1 1])<rntol)
        newburpar=1/6.*[-1 -2 1];
        newstairod=1/6.*[-1 1 0];
        splitdir=cross([-1 1 1]/norm([-1 1 1]),tangent)/norm(cross([-1 1 1]/norm([-1 1 1]),tangent));  
        
    elseif(norm(burgv-1/6.*[-2 -1 -1])<rntol)
        newburpar=1/6.*[-1 -2 -1];
        newstairod=1/6.*[-1 1 0];
        splitdir=cross([1 -1 1]/norm([1 -1 1]),tangent)/norm(cross([1 -1 1]/norm([1 -1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[2 1 1])<rntol)
        newburpar=1/6.*[1 2 1];
        newstairod=1/6.*[1 -1 0];
        splitdir=cross([1 -1 1]/norm([1 -1 1]),tangent)/norm(cross([1 -1 1]/norm([1 -1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[-1 -2 -1])<rntol)
        newburpar=1/6.*[-2 -1 -1];
        newstairod=1/6.*[1 -1 0];
        splitdir=cross([-1 1 1]/norm([-1 1 1]),tangent)/norm(cross([-1 1 1]/norm([-1 1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[1 2 1])<rntol)
        newburpar=1/6.*[2 1 1];
        newstairod=1/6.*[-1 1 0];
        splitdir=cross([-1 1 1]/norm([-1 1 1]),tangent)/norm(cross([-1 1 1]/norm([-1 1 1]),tangent)); 
        
    end
elseif((norm(cross(tangent/norm(tangent),[0 1 1]/norm([0 1 1])))<rntol)|(norm(cross(tangent/norm(tangent),[0 -1 -1]/norm([0 -1 -1])))<rntol))
    if(norm(burgv-1/6.*[-1 2 1])<rntol)
        newburpar=1/6.*[-1 1 2];
        newstairod=1/6.*[0 1 -1];
        splitdir=cross([1 -1 1]/norm([1 -1 1]),tangent)/norm(cross([1 -1 1]/norm([1 -1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[1 -2 -1])<rntol)
        newburpar=1/6.*[1 -1 -2];
        newstairod=1/6.*[0 -1 1];
        splitdir=cross([1 -1 1]/norm([1 -1 1]),tangent)/norm(cross([1 -1 1]/norm([1 -1 1]),tangent));   
        
    elseif(norm(burgv-1/6.*[-1 1 2])<rntol)
        newburpar=1/6.*[-1 2 1];
        newstairod=1/6.*[0 -1 1];
        splitdir=cross([1 1 -1]/norm([1 1 -1]),tangent)/norm(cross([1 1 -1]/norm([1 1 -1]),tangent));  
        
    elseif(norm(burgv-1/6.*[1 -1 -2])<rntol)
        newburpar=1/6.*[1 -2 -1];
        newstairod=1/6.*[0 1 -1];
        splitdir=cross([1 1 -1]/norm([1 1 -1]),tangent)/norm(cross([1 1 -1]/norm([1 1 -1]),tangent));  
        
    elseif(norm(burgv-1/6.*[-1 -1 -2])<rntol)
        newburpar=1/6.*[-1 -2 -1];
        newstairod=1/6.*[0 1 -1];
        splitdir=cross([1 -1 1]/norm([1 -1 1]),tangent)/norm(cross([1 -1 1]/norm([1 -1 1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[1 1 2])<rntol)
        newburpar=1/6.*[1 2 1];
        newstairod=1/6.*[0 -1 1];
        splitdir=cross([1 -1 1]/norm([1 -1 1]),tangent)/norm(cross([1 -1 1]/norm([1 -1 1]),tangent));  
        
    elseif(norm(burgv-1/6.*[-1 -2 -1])<rntol)
        newburpar=1/6.*[-1 -1 -2];
        newstairod=1/6.*[0 -1 1];
        splitdir=cross([1 1 -1]/norm([1 1 -1]),tangent)/norm(cross([1 1 -1]/norm([1 1 -1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[1 2 1])<rntol)
        newburpar=1/6.*[1 1 2];
        newstairod=1/6.*[0 1 -1];
        splitdir=cross([1 1 -1]/norm([1 1 -1]),tangent)/norm(cross([1 1 -1]/norm([1 1 -1]),tangent)); 
        
    end
    
elseif(norm(cross(tangent/norm(tangent),[1 0 1]/norm([1 0 1])))<rntol)  
    if(norm(burgv-1/6.*[2 -1 1])<rntol)
        newburpar=1/6.*[1 -1 2];
        newstairod=1/6.*[1 0 -1];
        splitdir=cross([1 1 -1]/norm([1 1 -1]),tangent)/norm(cross([1 1 -1]/norm([1 1 -1]),tangent)); 
        
    elseif(norm(burgv-1/6.*[-2 1 -1])<rntol)
        newburpar=1/6.*[-1 1 -2];
        newstairod=1/6.*[-1 0 1];
        splitdir=cross([1 1 -1]/norm([1 1 -1]),tangent)/norm(cross([1 1 -1]/norm([1 1 -1]),tangent));   
    elseif(norm(burgv-1/6.*[1 -1 2])<rntol)
        newburpar=1/6.*[2 -1 1];
        newstairod=1/6.*[-1 0 1];
        splitdir=cross([-1 1 1]/norm([-1 1 1]),tangent)/norm(cross([-1 1 1]/norm([-1 1 1]),tangent));   
    elseif(norm(burgv-1/6.*[-1 1 -2])<rntol)
        newburpar=1/6.*[-2 1 -1];
        newstairod=1/6.*[1 0 -1];
        splitdir=cross([-1 1 1]/norm([-1 1 1]),tangent)/norm(cross([-1 1 1]/norm([-1 1 1]),tangent));   
    elseif(norm(burgv-1/6.*[-1 -1 -2])<rntol)
        newburpar=1/6.*[-2 -1 -1];
        newstairod=1/6.*[1 0 -1];
        splitdir=cross([-1 1 1]/norm([-1 1 1]),tangent)/norm(cross([-1 1 1]/norm([-1 1 1]),tangent));   
    elseif(norm(burgv-1/6.*[1 1 2])<rntol)
        newburpar=1/6.*[2 1 1];
        newstairod=1/6.*[-1 0 1];
        splitdir=cross([-1 1 1]/norm([-1 1 1]),tangent)/norm(cross([-1 1 1]/norm([-1 1 1]),tangent));   
    elseif(norm(burgv-1/6.*[-2 -1 -1])<rntol)
        newburpar=1/6.*[-1 -1 -2];
        newstairod=1/6.*[-1 0 1];
        splitdir=cross([1 1 -1]/norm([1 1 -1]),tangent)/norm(cross([1 1 -1]/norm([1 1 -1]),tangent));   
    elseif(norm(burgv-1/6.*[2 1 1])<rntol)
        newburpar=1/6.*[1 1 2];
        newstairod=1/6.*[1 0 -1];
        splitdir=cross([1 1 -1]/norm([1 1 -1]),tangent)/norm(cross([1 1 -1]/norm([1 1 -1]),tangent));   
    end
end
    
      
     
     
     
     
     
     
     
     
     
     
    