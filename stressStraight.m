function [stress] = stressStraight(P0,P1, b, x)
%Function to return the stress at point x from a dislocation segment running 
%from vector b point p0 to p1 with burgers

length = norm(P0-P1); %length of dislocation segment
t = (P0-P1)/length; %normalized line tangent
bCt = cross(b,t); %plane normal

%Material Parameters
if ~exist('NU')
    NU=0.33;
end
C1 = 1-NU;
C2 = 1/(4*pi*C1);

tempStress = nonSymmStress(P0, P1, b, x); %Get the non-symmetrical stress

stress =  C2*(tempStress+transpose(tempStress));

end


function [nonSymKernel] = nonSymmStress_kernel(r,P0, P1, b)
%Cai's non-singular theory for dislocation stress

a2 = 2.32; %non-schmid coefficient for BCC metals
length = norm(P0-P1); %length of dislocation segment
t = (P0-P1)/length; %normalized line tangent
bCt = cross(b,t); %plane normal

%Material Parameters
if ~exist('NU')
    NU=0.33;
end

C1 = 1-NU;
C2 = 1/(4*pi*C1);

%Giacomo's code, for reference
%             const Scalar Ra2=r.squaredNorm()+DislocationStress<dim>::a2;
%             const Scalar Ra(sqrt(Ra2));
%             const VectorDim Ya(r+Ra*t);
%             const Scalar Yat(Ya.dot(t));
%             const Scalar Ya2a2(Ya.squaredNorm()+DislocationStress<dim>::a2);
%             const VectorDim bYa(b.cross(Ya));
%             const Scalar bYat(bYa.dot(t));
%                    
%             const Scalar f1(2.0/Ya2a2);

             Ra2=dot(r,r)+a2;
             Ra=sqrt(Ra2);
             Ya=r+Ra*t;
             Yat= dot(Ya,t);
             Ya2a2 = dot(Ya,Ya)+a2;
             bYa=cross(b,Ya);
             bYat=dot(bYa,t);
            
            
            f1=2.0/Ya2a2;
            
            nonSymKernel= f1*C1*(1.0+a2/Ya2a2)*t*transpose(bYa)+...
            f1*C1*0.5*a2/Ra2*t*transpose(cross(b,r))-...
            f1*Ya*transpose(bCt)-...
            f1*bYat/Yat*t*transpose(r)-...
            0.5*f1*bYat*(1.0+2.0*a2/Ya2a2+a2/Ra2)*[1 0 0;0 1 0;0 0 1]-...
            f1*bYat*(Ra+Yat)/Ra/Ya2a2*r*transpose(r)-0.5*f1*bYat*Ra/Yat*t*transpose(t);
end

function [nonSymStressMat] = nonSymmStress(P0, P1, b, x)
%Cai's non-singular theory for dislocation stress

    nonSymStressMat = nonSymmStress_kernel(P1-x,P0, P1, b)-nonSymmStress_kernel(P0-x,P0, P1, b);
end