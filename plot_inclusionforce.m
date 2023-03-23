%This function gives the force on the nodes of a segment due to the elastic inclusion depending
%on where is the segment

%p_file=fopen('Inclusion_force_screw','w');

%This function adds the force due to a elastic spherical inclusion

MU = 8.2588e-9;
NU = 0.305;%MD Mishin
inclusion_pos_rad=[0 0 0 30 8.2588e-9 0.305 -1e-9 1e-25 1e-25]; 
num_inclusions=size(inclusion_pos_rad,1);
YoungModulus_2=2*(1+NU)*MU;
B_2=-(YoungModulus_2*(1-NU))/((1+NU)*(1-(2*NU)));
a=(3/2)*((1/2)-NU);
b=(1-2*NU);
c=NU - (1/2);
e=1-NU;

%d=1-3*NU;

center_inclusion_i=inclusion_pos_rad(1:3);
R_i=inclusion_pos_rad(4);
MU_i=inclusion_pos_rad(5);
NU_i=inclusion_pos_rad(6);
rho_i=inclusion_pos_rad(7);
YoungModulus_i=2*(1+NU_i)*MU_i;
B_i=-(YoungModulus_i*(1-NU_i))/((1+NU_i)*(1-(2*NU_i)));
%alpha_i=6+2*NU_i;
%beta_i=5*((3/2) + 2*NU_i);
%gamma_i=1+5*NU_i;
%delta_i=5*((1/2) + NU_i);
alpha_i=1-2*NU_i;
beta_i=3*((1/4)-(1/2)*NU_i);
gamma_i=1+3*NU_i;
delta_i=((3/2)+3*NU_i);
b1=[0.49 0.51 0]/norm([0.49  0.51  0]);
x1_o=[10  0  0] + 1e-3*[1 1 0];
x2_o=[20 -10  0] + 1e-3*[1 1 0];
dir=[1 -1 0]/norm([1 -1 0]);
for i=1:30000
    x1=x1_o + i*5e-3*dir;
    x2=x2_o + i*5e-3*dir;
    
    zbase=(x2-x1)/norm((x2-x1));
    ybase=cross(zbase,x2)/norm(cross(zbase,x2));
    xbase=cross(ybase,zbase)/norm(cross(ybase,zbase));
    Base=[xbase',ybase',zbase'];
    b1_newbase=(inv(Base)*b1')';
    x1_newbase=(inv(Base)*x1')';
    x2_newbase=(inv(Base)*x2')';
    %distance_1(i)=x1_newbase(1);%norm(x1_newbase);
    %distance_2(i)=x2_newbase(1);%norm(x2_newbase);
    distance_1(i)=norm(x1_newbase);
    distance_2(i)=norm(x2_newbase);
    d=x1_newbase(3);
    H_j=x1_newbase(1);
    L_j=norm(x2-x1);
    norm_x1_newbase=norm(x1_newbase);
    norm_x2_newbase=norm(x2_newbase);
    if((norm(x1_newbase)<=R_i)&(norm(x2_newbase)<=R_i)) %The segment is completely inside of the inclusion
        [finclusion_1_i, finclusion_2_i]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
    elseif((norm(x1_newbase)>R_i)&(norm(x2_newbase)>R_i))
        if((R_i^2-x1_newbase(1)^2-x1_newbase(2)^2)<=0)
            [finclusion_1_i, finclusion_2_i]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
        elseif((R_i^2-x1_newbase(1)^2-x1_newbase(2)^2)>0)
            x_cut_1(1:2)=x1_newbase(1:2);
            x_cut_1(3)=sqrt(R_i^2-x1_newbase(1)^2-x1_newbase(2)^2);
            x_cut_2(1:2)=x1_newbase(1:2);
            x_cut_2(3)=-sqrt(R_i^2-x1_newbase(1)^2-x1_newbase(2)^2);
            if(((x_cut_1(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_1(3)>max(x1_newbase(3),x2_newbase(3))))&((x_cut_2(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_2(3)>max(x1_newbase(3),x2_newbase(3)))))
                [finclusion_1_i, finclusion_2_i]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
            else
                L_inside=norm(x_cut_2 - x_cut_1);
                L_outside=L_j - L_inside;
                [finclusion_1_in, finclusion_2_in]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                [finclusion_1_out, finclusion_2_out]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                finclusion_1_i=(L_inside/L_j)*finclusion_1_in + (L_outside/L_j)*finclusion_1_out;
                finclusion_2_i=(L_inside/L_j)*finclusion_2_in + (L_outside/L_j)*finclusion_2_out;
            end
        end
    else
        x_cut(1:2)=x1_newbase(1:2);
        x_cut(3)=sqrt(R_i^2-x1_newbase(1)^2-x1_newbase(2)^2);
        if((x_cut(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut(3)>max(x1_newbase(3),x2_newbase(3))))
            x_cut(3)=-sqrt(R_i^2-x1_newbase(1)^2-x1_newbase(2)^2);
        end
        if(norm(x1_newbase)<=R_i)
            L_inside=norm(x_cut-x1_newbase);
            L_outside=L_j - L_inside;
        elseif(norm(x2_newbase)<=R_i)
            L_inside=norm(x_cut-x2_newbase);
            L_outside=L_j - L_inside;
        end
        [finclusion_1_in, finclusion_2_in]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
        [finclusion_1_out, finclusion_2_out]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
        finclusion_1_i=(L_inside/L_j)*finclusion_1_in + (L_outside/L_j)*finclusion_1_out;
        finclusion_2_i=(L_inside/L_j)*finclusion_2_in + (L_outside/L_j)*finclusion_2_out;
    end
    
    finclusion_1=(Base*finclusion_1_i')';
    finclusion_2=(Base*finclusion_2_i')';
    
    %force_on_1(i)=norm(finclusion_1);
    %force_on_2(i)=norm(finclusion_2);
    
    %force_on_1(i)=dot(finclusion_1,[1 1 1]/norm([1 1 1]));
    %force_on_2(i)=dot(finclusion_2,[1 1 1]/norm([1 1 1]));
    
    force_on_1(i)=dot(finclusion_1,x1/norm(x1));
    force_on_2(i)=dot(finclusion_2,x2/norm(x2));
    
    %force_on_1(i)=dot(finclusion_1,dir);
    %force_on_2(i)=dot(finclusion_2,dir);
    
    %fprintf(p_file,'%f\t%e\t%f\t%e\n',distance_1(i),force_on_1(i),distance_2(i),force_on_2(i));

end
%semilogx(distance_1,force_on_1,'r-',distance_2,force_on_2,'b-');
plot(distance_1,force_on_1,'r-',distance_1,force_on_2,'b-');
