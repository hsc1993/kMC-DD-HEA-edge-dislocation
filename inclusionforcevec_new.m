% This function gives the force on the nodes of a segment due to the elastic inclusion depending
% on where is the segment


function [finclusion_1, finclusion_2]=inclusionforcevec_new(MU,NU,segments,linkid,inclusion_pos_rad)

% This function adds the force due to a elastic spherical inclusion

lseg=size(segments,1);
num_inclusions=size(inclusion_pos_rad,1);
YoungModulus_2=2*(1+NU)*MU;
B_2=-(YoungModulus_2*(1-NU))/((1+NU)*(1-(2*NU)));


% d=1-3*NU;
finclusion_1=zeros(lseg,3);
finclusion_2=zeros(lseg,3);
if length(linkid)==0
    
    for i=1:num_inclusions
        finclusion_1_i=zeros(lseg,3);
        finclusion_2_i=zeros(lseg,3);
        center_inclusion_i=inclusion_pos_rad(i,1:3);
        R_i=inclusion_pos_rad(i,4);
        MU_i=inclusion_pos_rad(i,5);
        NU_i=inclusion_pos_rad(i,6);
        beta=inclusion_pos_rad(i,7); % Fitting parameter
        YoungModulus_i=2*(1+NU_i)*MU_i;
        B_i=-(YoungModulus_i*(1-NU_i))/((1+NU_i)*(1-(2*NU_i)));
        c=R_i*0.1; % The strain value at particle interface. Needs to be a function of R_i (and/or beta)!!!!! TODO!
        Rs=R_i*15;  % The radius where the strain is zero. Needs to be a function of R_i (and/or beta)!!!!! TODO!
        for j=1:lseg
            b1=segments(j,3:5);
            x1=segments(j,6:8)-inclusion_pos_rad(i,1:3);
            x2=segments(j,9:11)-inclusion_pos_rad(i,1:3);
            zbase=(x2-x1)/norm((x2-x1));
            ybase=cross (zbase,x2)/norm(cross(zbase,x2));
            xbase=cross (ybase,zbase)/norm(cross(ybase,zbase));
            Base=[xbase',ybase',zbase'];
            b1_newbase=(inv(Base)*b1')';
            x1_newbase=(inv(Base)*x1')';
            x2_newbase=(inv(Base)*x2')';
            d=x1_newbase(3);
            H_j=x1_newbase(1);
            L_j=norm(x2-x1);
            norm_x1_newbase=norm(x1_newbase);
            norm_x2_newbase=norm(x2_newbase);
            if((norm(x1_newbase)<=R_i)&(norm(x2_newbase)<=R_i)) % The segment is completely inside of the inclusion
                [finclusion_1_i(j,:), finclusion_2_i(j,:)]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta);
            elseif((norm(x1_newbase)>R_i)&(norm(x2_newbase)>R_i))
                if((R_i^2-x1_newbase (1)^2-x1_newbase (2)^2)<=0)
                    [finclusion_1_i(j,:), finclusion_2_i(j,:)]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta);
                elseif((R_i^2-x1_newbase (1)^2-x1_newbase (2)^2)>0)
                    x_cut_1(1:2)=x1_newbase(1:2);
                    x_cut_1(3)=sqrt(R_i^2-x1_newbase (1)^2-x1_newbase (2)^2);
                    x_cut_2(1:2)=x1_newbase(1:2);
                    x_cut_2(3)=-sqrt(R_i^2-x1_newbase (1)^2-x1_newbase (2)^2);
                    if(((x_cut_1(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_1(3)>max(x1_newbase(3),x2_newbase(3))))&((x_cut_2(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_2(3)>max(x1_newbase(3),x2_newbase(3)))))
                        [finclusion_1_i(j,:), finclusion_2_i(j,:)]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta);
                    else
                        L_inside=norm(x_cut_2 - x_cut_1);
                        L_outside=L_j - L_inside;
                        [finclusion_1_in, finclusion_2_in]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta);
                        [finclusion_1_out, finclusion_2_out]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta);
                        finclusion_1_i(j,:)=(L_inside/L_j)*finclusion_1_in + (L_outside/L_j)*finclusion_1_out;
                        finclusion_2_i(j,:)=(L_inside/L_j)*finclusion_2_in + (L_outside/L_j)*finclusion_2_out;
                    end
                end
            else
                x_cut(1:2)=x1_newbase(1:2);
                x_cut(3)=sqrt(R_i^2-x1_newbase (1)^2-x1_newbase (2)^2);
                if((x_cut(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut(3)>max(x1_newbase(3),x2_newbase(3))))
                    x_cut(3)=-sqrt(R_i^2-x1_newbase (1)^2-x1_newbase (2)^2);
                end
                if(norm(x1_newbase)<=R_i)
                    L_inside=norm(x_cut-x1_newbase);
                    L_outside=L_j - L_inside;
                elseif(norm(x2_newbase)<=R_i)
                    L_inside=norm(x_cut-x2_newbase);
                    L_outside=L_j - L_inside;
                end
                [finclusion_1_in, finclusion_2_in]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta);
                [finclusion_1_out, finclusion_2_out]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta);
                finclusion_1_i(j,:)=(L_inside/L_j)*finclusion_1_in + (L_outside/L_j)*finclusion_1_out;
                finclusion_2_i(j,:)=(L_inside/L_j)*finclusion_2_in + (L_outside/L_j)*finclusion_2_out;
            end
            
        end
        finclusion_1=finclusion_1 + finclusion_1_i;
        finclusion_2=finclusion_2 + finclusion_2_i;
    end
    finclusion_1=(Base*finclusion_1')';
    finclusion_2=(Base*finclusion_2')';

else
    for i=1:num_inclusions
        finclusion_1_i=zeros(lseg,3);
        finclusion_2_i=zeros(lseg,3);
        center_inclusion_i=inclusion_pos_rad(i,1:3);
        R_i=inclusion_pos_rad(i,4);
        MU_i=inclusion_pos_rad(i,5);
        NU_i=inclusion_pos_rad(i,6);
        rho_i=inclusion_pos_rad(i,7);
        YoungModulus_i=2*(1+NU_i)*MU_i;
        B_i=-(YoungModulus_i*(1-NU_i))/((1+NU_i)*(1-(2*NU_i)));
        alpha_i=1-2*NU_i;
        beta_i=3*((1/4) - (1/2)*NU_i);
        gamma_i=1+3*NU_i;
        delta_i=(3/2) + 3*NU_i;
        for j=1:length(linkid)
            b1=segments(linkid(j),3:5);
            x1=segments(linkid(j),6:8)-inclusion_pos_rad(i,1:3);
            x2=segments(linkid(j),9:11)-inclusion_pos_rad(i,1:3);
            zbase=(x2-x1)/norm((x2-x1));
            ybase=cross (zbase,x2)/norm(cross(zbase,x2));
            xbase=cross (ybase,zbase)/norm(cross(ybase,zbase));
            Base=[xbase',ybase',zbase'];
            b1_newbase=(inv(Base)*b1')';
            x1_newbase=(inv(Base)*x1')';
            x2_newbase=(inv(Base)*x2')';
            d=x1_newbase(3);
            H_j=x1_newbase(1);
            L_j=norm(x2-x1);
            norm_x1_newbase=norm(x1_newbase);
            norm_x2_newbase=norm(x2_newbase);
            if((norm(x1_newbase)<=R_i)&(norm(x2_newbase)<=R_i)) % The segment is completely inside of the inclusion
                [finclusion_1_i(linkid(j),:), finclusion_2_i(linkid(j),:)]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
            elseif((norm(x1_newbase)>R_i)&(norm(x2_newbase)>R_i))
                if((R_i^2-x1_newbase (1)^2-x1_newbase (2)^2)<=0)
                    [finclusion_1_i(linkid(j),:), finclusion_2_i(linkid(j),:)]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                elseif((R_i^2-x1_newbase (1)^2-x1_newbase (2)^2)>0)
                    x_cut_1(1:2)=x1_newbase(1:2);
                    x_cut_1(3)=sqrt(R_i^2-x1_newbase (1)^2-x1_newbase (2)^2);
                    x_cut_2(1:2)=x1_newbase(1:2);
                    x_cut_2(3)=-sqrt(R_i^2-x1_newbase (1)^2-x1_newbase (2)^2);
                    if(((x_cut_1(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_1(3)>max(x1_newbase(3),x2_newbase(3))))&((x_cut_2(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_2(3)>max(x1_newbase(3),x2_newbase(3)))))
                        [finclusion_1_i(linkid(j),:), finclusion_2_i(linkid(j),:)]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                    else
                        L_inside=norm(x_cut_2 - x_cut_1);
                        L_outside=L_j - L_inside;
                        [finclusion_1_in, finclusion_2_in]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,alpha_i,beta_i,gamma_i,delta_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                        [finclusion_1_out, finclusion_2_out]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,a,b,c,e,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j);
                        finclusion_1_i(linkid(j),:)=(L_inside/L_j)*finclusion_1_in + (L_outside/L_j)*finclusion_1_out;
                        finclusion_2_i(linkid(j),:)=(L_inside/L_j)*finclusion_2_in + (L_outside/L_j)*finclusion_2_out;
                    end
                end
            else
                x_cut(1:2)=x1_newbase(1:2);
                x_cut(3)=sqrt(R_i^2-x1_newbase (1)^2-x1_newbase (2)^2);
                if((x_cut(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut(3)>max(x1_newbase(3),x2_newbase(3))))
                    x_cut(3)=-sqrt(R_i^2-x1_newbase (1)^2-x1_newbase (2)^2);
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
                finclusion_1_i(linkid(j),:)=(L_inside/L_j)*finclusion_1_in + (L_outside/L_j)*finclusion_1_out;
                finclusion_2_i(linkid(j),:)=(L_inside/L_j)*finclusion_2_in + (L_outside/L_j)*finclusion_2_out;
            end
            
        end
        finclusion_1=finclusion_1 + finclusion_1_i;
        finclusion_2=finclusion_2 + finclusion_2_i;
    end
    finclusion_1=(Base*finclusion_1')';
    finclusion_2=(Base*finclusion_2')';
    
end

function [finclusion_1, finclusion_2]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta)

    b1=b1_newbase(1);
    b2=b1_newbase(2);
    b3=b1_newbase(3);
    
    x=x1_newbase(1);
    z_1=x1_newbase(3);
    z_2=x2_newbase(3);
    
	B=B_i;
	aI=(NU_i-3);
	bI=(1+NU_i);
	cI=(2*NU_i-1);
	dI=(1+2*NU_i);
	eI=(2+6*NU_i);
	fI=(2*NU_i-3);
	G=c*beta/(R_i-Rs)*exp((R_i-Rs)/beta);
	rhoI=(10*G*R_i-20*c)/(B*R_i^2);
    const_general=MU_i/(120*R_i*bI*cI);
    
    
        %%%%%%%%%%
        %
        %   This section calculates the total force on the segment due to the inclusion
        %
        %%%%%%%%%%

        
    finclusion_total=zeros(1,3);
    
    finclusion_total(1) = const_general*(-12*aI*B*b2*H_j^2*rhoI*z1 - ...
	12*B*b3*cI*H_j*rhoI*z_1^2 + 2*B*b2*eI*rhoI*z1^3 - ...
	3*(10*G - 3*B*R_i*rho_I)*(2*b3*cI*H_j - b2*dI*z_1)*sqrt(H_j^2 + z_1^2) + ...
	12*aI*B*b2*H_j^2*rhoI*z_2 + 12*B*b3*cI*H_j*rhoI*z_2^2 - 2*B*b2*eI*rhoI*z_2^3 + ...
	3*(10*G - 3*B*R_i*rhoI)*(2*b3*cI*H_j - b2*dI*z_2)*sqrt(H_j^2 + z_2^2) - ...
	3*b2*fI*H_j^2*(10*G - 3*B*R_i*rhoI)*log(z_1 + sqrt(H_j^2 + z_1^2)) + ...
	3*b2*fI*H_j^2*(10*G - 3*B*R_i*rhoI)*log(z_2 + sqrt(H_j^2 + z_2^2)));
    
    finclusion_total(2) = const_general*b1*(-6*B*eI*H_j^2*rhoI*z_1 - 2*B*eI*rhoI*z_1^3 - ...
	30*dI*G*z_1*sqrt(H_j^2 + z_1^2) + 9*B*dI*R_i*rhoI*z_1*sqrt(H_j^2 + z_1^2) + ...
	6*B*eI*H^2*rhoI*z_2 + 2*B*eI*rhoI*z_2^3 + ... 
	30*dI*G*z_2*sqrt(H_j^2 + z_2^2) - 9*B*dI*R_i*rhoI*z_2*sqrt(H^2 + z_2^2) + ...
	3*dI*H_j^2*(-10*G + 3*B*R_i*rhoI)*log(z_1 + sqrt(H_j^2 + z_1^2)) + ...
	3*dI*H_j^2*(10*G - 3*B*R_i*rhoI)*log(z_2 + sqrt(H_j^2 + z_2^2)));

    finclusion_total(3) = 0;

    %%%%%%%%%%
    %
    %   This section calculates the force on node 2
    %
    %%%%%%%%%%


    finclusion_2=zeros(1,3);
    
    finclusion_2(1) = const_general*((1/(2*L_j))*(-12*aI*B*b2*H_j^2*rhoI*z_1^2 - ...
	16*B*b3*cI*H_j*rhoI*z_1^3 + 3*B*b2*eI*rhoI*z_1^4 - ...
	2*(10*G - 3*B*R_i*rhoI)*sqrt(H_j^2 + z_1^2)*(3*b3*cI*H_j*z_1 + 2*b2*(2*(-3 + dI)*H_j^2 - dI*z_1^2)) + ...
	12*aI*B*b2*H_j^2*rhoI*z_2^2 + 16*B*b3*cI*H_j*rhoI*z2^3 - 3*B*b2*eI*rhoI*z_2^4 + ... 
	2*(10*G - 3*B*R_i*rhoI)*sqrt(H_j^2 + z_2^2)*(3*b3*cI*H_j*z_2 + 2*b2*(2*(-3 + dI)*H_j^2 - dI*z_2^2)) + ...
	6*b3*cI*H_j^3*(10*G - 3*B*R_i*rhoI)*log(z_1 + sqrt(H_j^2 + z_1^2)) - ...
	6*b3*cI*H_j^3*(10*G - 3*B*R_i*rhoI)*log(z_2 + sqrt(H_j^2 + z_2^2))));
    
    finclusion_2(2) = const_general*(b1/(2*L_j)*(-4*dI*(10*G - 3*B*R_i*rhoI)*(H_j^2 + z_1^2)^(3/2) - ...
	3*B*eI*rhoI*(H_j^2 + z_1^2)^2 + 4*dI*(10*G - 3*B*R_i*rhoI)*(H_j^2 + z_2^2)^(3/2) + 3*B*eI*rhoI*(H_j^2 + z_2^2)^2));

    finclusion_2(3) = 0;


    %%%%%%%%%%
    %
    %   This section calculates the force on node 1
    %
    %%%%%%%%%%

    finclusion_1=finclusion_total - finclusion_2;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [finclusion_1, finclusion_2]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,rho_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta,YoungModulus_2)


    
    b1=b1_newbase(1);
    b2=b1_newbase(2);
    b3=b1_newbase(3);
    
    x=x1_newbase(1);
    z_1=x1_newbase(3);
    z_2=x2_newbase(3);

    
    const_general=c*exp(Rs/beta)*YoungModulus_2/((-R + Rs)*(1 + NU));

    %%%%%%%%%%
    %
    %   This section calculates the total force on the segment due to the inclusion
    %
    %%%%%%%%%%

    
    finclusion_total=zeros(1,3);
    
    f_total_1(1)= b2*H_j^2/(2*beta^2)*(-((exp(-(sqrt (H^2 + z1^2)/beta))*(beta*R_i*(beta + sqrt(H^2 + z_1^2)) + ...
	exp(sqrt (H_j^2 + z_1^2)/beta)*(2*beta + R_i)*(H_j^2 + z_1^2)*expint(-(sqrt (H_j^2 + z_1^2)/beta))))/(H_j^2 + z_1^2)) + ...
	(exp(-(sqrt (H_j^2 + z_2^2)/beta))*(beta*R_i*(beta + sqrt(H_j^2 + z_2^2)) + ...
	exp(sqrt (H_j^2 + z_2^2)/beta)*(2*beta + R_i)*(H_j^2 + z_2^2)*expint(-(sqrt (H_j^2 + z_2^2)/beta))))/(H_j^2 + z_2^2));

    f_total_2(1)=b1/(-1 + 2*NU)*(exp(-(sqrt (H_j^2 + z_1^2)/beta))*(R_i + beta*NU - sqrt(H_j^2 + z_1^2) + ...
	exp(sqrt (H_j^2 + z_1^2)/beta)*R_i*NU*expint(-(sqrt (H_j^2 + z_1^2)/beta))) - ...
	exp(-(sqrt (H_j^2 + z_2^2)/beta))*(R_i + beta*NU - sqrt(H_j^2 + z_2^2) + ...
	exp(sqrt (H_j^2 + z_2^2)/beta)*R_i*NU*expint(-(sqrt (H_j^2 + z_2^2)/beta))));

    f_total_3(1)=b3*H_j*((exp(-(sqrt (H_j^2 + z_1^2)/beta))*sqrt(z_1^2)*(R_i - sqrt(H_j^2 + z_1^2)))/(z_1*sqrt(H_j^2 + z_1^2)) + ...
	(exp(-(sqrt (H_j^2 + z_2^2)/beta))*sqrt(z_2^2)*(-R_i + sqrt(H_j^2 + z_2^2)))/(z_2*sqrt(H^2 + z2^2)));
    
    finclusion_total(1)=const_general*(f_total_1(1) + f_total_2(1) + f_total_3(1));
    
    finclusion_total(2)=-const_general*(b1/(-1 + 2*NU)*(exp(-(sqrt(H_j^2 + z_1^2)/beta))*(R_i + beta*NU - sqrt(H_j^2 + z_1^2) + ...
	(exp((sqrt(H_j^2 + z_1^2)/beta))*R_i*NU*expint(-(sqrt(H_j^2 + z_1^2)/beta))) - ...
	(exp(-(sqrt(H_j^2 + z_2^2)/beta))*(R_i + beta*NU - sqrt(H_j^2 + z_2^2) + ...
	(exp((sqrt(H_j^2 + z_2^2)/beta))*R_i*NU*expint(-(sqrt(H_j^2 + z_2^2)/beta))))))));
	
    finclusion_total(3)=0;

    %%%%%%%%%%
    %
    %   This section calculates the force on node 2
    %
    %%%%%%%%%%
    
    
    finclusion_2=zeros(1,3);
    
    f_2_1(1) = ((b2*H_j^2)/(L_j)*(exp(-(sqrt(H_j^2 + z_1^2)/beta))*(-1 + R_i/sqrt(H_j^2 + z_1^2)) + ...
	exp(-(sqrt(H_j^2 + z_2^2)/beta))*(1 - R_i/sqrt(H_j^2 + z_2^2))));
	
    f_2_2(1) = (b2/(L_j*(-1 + 2*NU))*(exp(-(sqrt(H_j^2 + z_1^2)/beta))*(H_j^2 - beta*(beta - R_i)*(-1 + NU) + ...
	z_1^2 - (R_i + beta*(-1 + NU))*sqrt(H_j^2 + z_1^2)) + ...
	exp(-(sqrt(H_j^2 + z_2^2)/beta))*(-H_j^2 + beta*(beta - R_i)*(-1 + NU) - ...
	z_2^2 + (R_i + beta*(-1 + NU))*sqrt(H_j^2 + z_2^2))));
	
    f_2_3(1)=(b3 *H_j/(2*beta^2*L_j)*(-((exp(-(sqrt(H_j^2 + z_1^2)/beta))*(-beta (H_j^2*R_i*sqrt(H_j^2 + z_1^2) +...
	2*beta^2*(H_j^2 + z_1^2) + beta*(-H_j^2*R_i - 2*R_i*z_1^2 + 2*(H_j^2 + z_1^2)^(3/2))) - ...
	exp(sqrt(H_j^2 + z_1^2)/beta)*(2*beta*H_j^2 + 2*beta^2*R_i + H_j^2*R_i)*(H_j^2 + z_1^2)*expint(-(sqrt(H_j^2 + z_1^2)/beta))))/(H_j^2 + z_1^2)) + ...
	(exp(-(sqrt(H_j^2 + z_2^2)/beta))*(-beta*(H_j^2*R_i*sqrt(H_j^2 + z_2^2) + 2*beta^2*(H_j^2 + z_2^2) + beta*(-H_j^2*R - 2*R*z_2^2 + 2*(H_j^2 + z_2^2)^(3/2))) - ... 
	exp(sqrt(H_j^2 + z_2^2)/beta)*(2*beta*H_j^2 + 2*beta^2*R_i + H_j^2*R)*(H_j^2 + z_2^2)*expint(-(sqrt(H_j^2 + z_2^2)/beta))))/(H_j^2 + z_2^2)));

    finclusion_2(1)=const_general*(f_2_1(1) + f_2_2(1) + f_2_3(1));
    
    finclusion_2(2)=-const_general*(b1/(L_j*(-1 + 2*NU))*(exp(-(sqrt(H_j^2 + z_1^2)/beta))*(-H_j^2 + beta*(beta - R_i)*(-1 + NU) - ...
	z_1^2 + (R_i + beta*(-1 + NU))*sqrt(H_j^2 + z_1^2)) + exp(-(sqrt(H_j^2 + z_2^2)/beta)) (H_j^2 - beta*(beta - R_i)*(-1 + NU) + ... 
	z_2^2 - (R_i + beta (-1 + NU))*sqrt(H_j^2 + z_2^2))));
    
    finclusion_2(3)=0;

    %%%%%%%%%%
    %
    %   This section calculates the force on node 1
    %
    %%%%%%%%%%

    finclusion_1=finclusion_total - finclusion_2;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
