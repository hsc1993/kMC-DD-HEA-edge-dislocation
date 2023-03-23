%This function gives the force on the nodes of a segment due to the elastic inclusion depending
%on where is the segment


function [finclusion_1, finclusion_2]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j)


    eps=1e-6;
    
    finclusion_1=zeros(1,3);
    finclusion_2=zeros(1,3);
    
    b1=b1_newbase(1);
    b2=b1_newbase(2);
    b3=b1_newbase(3);
    
    x=x1_newbase(1);
    z_1=x1_newbase(3);
    z_2=x2_newbase(3);
    
    norm_r1=sqrt(x^2+z_1^2);
    norm_r2=sqrt(x^2+z_2^2);
    
    rho_o=(1/5)*rho_i*B_i/B_2;
    %rho_o=(3/14)*rho_i*B_i/B_2;
    
    const_general=rho_o*R_i^3/(1-NU);

    %%%%%%%%%%
    %
    %   This section calculates the total force on the segment due to the inclusion
    %
    %%%%%%%%%%

    integral_node_1=zeros(1,3);
    integral_node_2=zeros(1,3);
    finclusion_total=zeros(1,3);
    
    f_total_1(1)=b2*H_j^2*((R_i*alpha_i/(3*x^4))*((2*z_2^3 + 3*x^2*z_2)/norm_r2^3 - (2*z_1^3 + 3*x^2*z_1)/norm_r1^3) - (beta_i/(2*x^3))*(x*(z_2/norm_r2^2 - z_1/norm_r1^2) + (atan(z_2/x) - atan(z_1/x))));
    f_total_2(1)=(b2/2)*((R_i*gamma_i/x^2)*(z_2/norm_r2 - z_1/norm_r1) + (delta_i/x)*(atan(z_2/x) - atan(z_1/x)));
    f_total_3(1)=b3*H_j*(-(R_i*alpha_i/3)*((1/norm_r2^(3)) - (1/norm_r1^(3))) + (beta_i/2)*((1/norm_r2^2) - (1/norm_r1^2)));
    
    finclusion_total(1)=const_general*(f_total_1(1) + f_total_2(1) + f_total_3(1));
    
    f_total_1(2)=-(b1/2)*((R_i*gamma_i/x^2)*(z_2/norm_r2 - z_1/norm_r1) + (delta_i/x)*(atan(z_2/x) - atan(z_1/x)));
    finclusion_total(2)=const_general*f_total_1(2);
    
    finclusion_total(3)=0;

    %%%%%%%%%%
    %
    %   This section calculates the force on node 2
    %
    %%%%%%%%%%
    
    %f_2_1(1)=-(d/L_j)*f_total_1(1) + b2*(H_j^2/L_j)*(-(R_i*alpha_i/3)*(1/norm_r2^3 - 1/norm_r1^3) + (beta_i/2)*(1/norm_r2^2 - 1/norm_r1^2));
    %f_2_2(1)=-(d/L_j)*f_total_2(1) + (b2/2)*(1/L_j)*(-R_i*gamma_i*(1/norm_r2 - 1/norm_r1) + (delta_i/2)*log(norm_r2^2/norm_r1^2));  
    %f_2_3(1)=-(d/L_j)*f_total_3(1) + b3*(H_j/L_j)*(R_i*alpha_i/(3*x^2)*(z_2^3/norm_r2^(3) - z_1^3/norm_r1^(3)) -(beta_i/2)*((atan(z_2/x) - atan(z_1/x))/x - (z_2/norm_r2^2 - z_1/norm_r1^2)));
    
    %finclusion_2(1)=const_general*(f_2_1(1) + f_2_2(1) + f_2_3(1));
    
    %f_2_1(2)=-(d/L_j)*f_total_1(2) - (b1/2)*(1/L_j)*(-R_i*gamma_i*(1/norm_r2 - 1/norm_r1) + (delta_i/2)*log(norm_r2^2/norm_r1^2));
    %finclusion_2(2)=const_general*(b1/2)*f_2_1(2);
    
    %finclusion_2(3)=0;
    finclusion_2=finclusion_total/2;
    %%%%%%%%%%
    %
    %   This section calculates the force on node 1
    %
    %%%%%%%%%%

    %finclusion_1=finclusion_total - finclusion_2;
    finclusion_1=finclusion_total/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
