function [finclusion_1, finclusion_2]=test_inclusionforcevec_new()

steps = 60;
nf1 = zeros(steps,1);
nf2 = zeros(steps,1);
disp = zeros(steps,1);
step = 0.5;
f1=zeros(steps,3);
f2=zeros(steps,3);

for i = 1:steps
    rn    = [
            2.33+(i-1)*step   1+(i-1)*step   0 7;
             2.33+(i-1)*step   -1+(i-1)*step  0 7;      
             
         ];
%     plot(rn(:,1),rn(:,2));
%     hold on;
%     norm(rn(1,1:3)-rn(2,1:3))
    
    b1 = [  1 -1 0 ]/2;
    n1 = [  -1 -1 -1 ]/sqrt(3);

    links = [
             1   2  b1 n1;
           ];

    linkid = [];

    segments=constructsegmentlist(rn,links);   

    MU=4.5628*10^-9;
    NU=1;
    Biscrew=1;
    Biedge=1;
    MUi=5.4461*10^-9;
    NUi=0.274;
    betai=28.35874;
    radius=2.33;
    inclusion_pos_rad=[0 0 0 radius MUi NUi betai Biscrew Biedge];

    [finclusion_1, finclusion_2]=inclusionforcevec_new(MU,NU,segments,linkid,inclusion_pos_rad);
    f1(i,:) = finclusion_1;
    f2(i,:) = finclusion_2;
    nf1(i) = norm(finclusion_1);
    nf2(i) = norm(finclusion_2);
    disp(i) = norm(rn(1,1:3)-inclusion_pos_rad(1:3));
end
size(disp)
plot(disp,f1(:,1));
hold on
plot(disp,f1(:,2));
plot(disp,f1(:,3));
plot(disp,f2(:,1));
plot(disp,f2(:,2));
plot(disp,f2(:,3));

% theta = linspace(0,2*pi);
% x = radius*cos(theta);
% y = radius*sin(theta) ;
% plot(x,y)
hold off
legend('f node 1 in x','f node 1 in y','f node 1 in z','f node 2 in x','f node 2 in y','f node 2 in z');


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
        c=0.17; % The strain value at particle interface. Needs to be a function of R_i (and/or beta)!!!!! TODO!
        Rs=30.69265;  % The radius where the strain is zero. Needs to be a function of R_i (and/or beta)!!!!! TODO!
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
            
            %disp(sprintf('H_j=%d \n',H_j));
            %disp(sprintf('L_j=%d  \n',L_j));
            %disp(sprintf('d=%d  \n',d));
            norm_x1_newbase=norm(x1_newbase);
            norm_x2_newbase=norm(x2_newbase);
            if((norm(x1_newbase)<=R_i)&(norm(x2_newbase)<=R_i)) % The segment is completely inside of the inclusion
                [finclusion_1_i(j,:), finclusion_2_i(j,:)]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta);
            elseif((norm(x1_newbase)>R_i)&(norm(x2_newbase)>R_i))
                if((R_i^2-x1_newbase (1)^2-x1_newbase (2)^2)<=0)  % The segment is completely outside of the inclusion
                    [finclusion_1_i(j,:), finclusion_2_i(j,:)]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta,YoungModulus_2);
                elseif((R_i^2-x1_newbase (1)^2-x1_newbase (2)^2)>0)
                    x_cut_1(1:2)=x1_newbase(1:2);
                    x_cut_1(3)=sqrt(R_i^2-x1_newbase (1)^2-x1_newbase (2)^2);
                    x_cut_2(1:2)=x1_newbase(1:2);
                    x_cut_2(3)=-sqrt(R_i^2-x1_newbase (1)^2-x1_newbase (2)^2);
                    if(((x_cut_1(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_1(3)>max(x1_newbase(3),x2_newbase(3))))&((x_cut_2(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_2(3)>max(x1_newbase(3),x2_newbase(3)))))
                        [finclusion_1_i(j,:), finclusion_2_i(j,:)]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta,YoungModulus_2);
                    else % The segment is both outside and inside of the inclusion
                        L_inside=norm(x_cut_2 - x_cut_1);
                        L_outside=L_j - L_inside;
                        [finclusion_1_in, finclusion_2_in]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta);
                        [finclusion_1_out, finclusion_2_out]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta,YoungModulus_2);
                        finclusion_1_i(j,:)=(L_inside/L_j)*finclusion_1_in + (L_outside/L_j)*finclusion_1_out;
                        finclusion_2_i(j,:)=(L_inside/L_j)*finclusion_2_in + (L_outside/L_j)*finclusion_2_out;
                    end
                end
            else % The segment is both outside and inside of the inclusion
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
                else
                    L_inside=norm(x_cut-x2_newbase);
                    L_outside=L_j - L_inside;
                end
                [finclusion_1_in, finclusion_2_in]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta);
                [finclusion_1_out, finclusion_2_out]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta,YoungModulus_2);
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
        rho_i=inclusion_pos_rad(i,7); % Fitting parameter
        beta=inclusion_pos_rad(i,7); % Fitting parameter
        YoungModulus_i=2*(1+NU_i)*MU_i;
        B_i=-(YoungModulus_i*(1-NU_i))/((1+NU_i)*(1-(2*NU_i)));
        c=R_i*0.1; % The strain value at particle interface. Needs to be a function of R_i (and/or beta)!!!!! TODO!
        Rs=R_i*15;  % The radius where the strain is zero. Needs to be a function of R_i (and/or beta)!!!!! TODO!
        
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
                [finclusion_1_i(linkid(j),:), finclusion_2_i(linkid(j),:)]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta);
            elseif((norm(x1_newbase)>R_i)&(norm(x2_newbase)>R_i))
                if((R_i^2-x1_newbase (1)^2-x1_newbase (2)^2)<=0)
                    [finclusion_1_i(linkid(j),:), finclusion_2_i(linkid(j),:)]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta,YoungModulus_2);
                elseif((R_i^2-x1_newbase (1)^2-x1_newbase (2)^2)>0)
                    x_cut_1(1:2)=x1_newbase(1:2);
                    x_cut_1(3)=sqrt(R_i^2-x1_newbase (1)^2-x1_newbase (2)^2);
                    x_cut_2(1:2)=x1_newbase(1:2);
                    x_cut_2(3)=-sqrt(R_i^2-x1_newbase (1)^2-x1_newbase (2)^2);
                    if(((x_cut_1(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_1(3)>max(x1_newbase(3),x2_newbase(3))))&((x_cut_2(3)<min(x1_newbase(3),x2_newbase(3)))|(x_cut_2(3)>max(x1_newbase(3),x2_newbase(3)))))
                        [finclusion_1_i(linkid(j),:), finclusion_2_i(linkid(j),:)]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta,YoungModulus_2);
                    else
                        L_inside=norm(x_cut_2 - x_cut_1);
                        L_outside=L_j - L_inside;
                        [finclusion_1_in, finclusion_2_in]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta);
                        [finclusion_1_out, finclusion_2_out]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta,YoungModulus_2);
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
                [finclusion_1_in, finclusion_2_in]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta);
                [finclusion_1_out, finclusion_2_out]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta,YoungModulus_2);
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

%This function gives the force on the nodes of a segment due to the elastic inclusion depending
%on where is the segment

function [finclusion_1, finclusion_2]=ForceInsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta)

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
    
    finclusion_total(1) = const_general*(-12*aI*B*b2*H_j^2*rhoI*z_1 - ...
	12*B*b3*cI*H_j*rhoI*z_1^2 + 2*B*b2*eI*rhoI*z_1^3 - ...
	3*(10*G - 3*B*R_i*rhoI)*(2*b3*cI*H_j - b2*dI*z_1)*sqrt(H_j^2 + z_1^2) + ...
	12*aI*B*b2*H_j^2*rhoI*z_2 + 12*B*b3*cI*H_j*rhoI*z_2^2 - 2*B*b2*eI*rhoI*z_2^3 + ...
	3*(10*G - 3*B*R_i*rhoI)*(2*b3*cI*H_j - b2*dI*z_2)*sqrt(H_j^2 + z_2^2) - ...
	3*b2*fI*H_j^2*(10*G - 3*B*R_i*rhoI)*log(z_1 + sqrt(H_j^2 + z_1^2)) + ...
	3*b2*fI*H_j^2*(10*G - 3*B*R_i*rhoI)*log(z_2 + sqrt(H_j^2 + z_2^2)));
    
    finclusion_total(2) = const_general*b1*(-6*B*eI*H_j^2*rhoI*z_1 - 2*B*eI*rhoI*z_1^3 - ...
	30*dI*G*z_1*sqrt(H_j^2 + z_1^2) + 9*B*dI*R_i*rhoI*z_1*sqrt(H_j^2 + z_1^2) + ...
	6*B*eI*H_j^2*rhoI*z_2 + 2*B*eI*rhoI*z_2^3 + ... 
	30*dI*G*z_2*sqrt(H_j^2 + z_2^2) - 9*B*dI*R_i*rhoI*z_2*sqrt(H_j^2 + z_2^2) + ...
	3*dI*H_j^2*(-10*G + 3*B*R_i*rhoI)*log(z_1 + sqrt(H_j^2 + z_1^2)) + ...
	3*dI*H_j^2*(10*G - 3*B*R_i*rhoI)*log(z_2 + sqrt(H_j^2 + z_2^2)));

    finclusion_total(3) = 0;
    
    %disp(sprintf('finclusion_total outside=%d %d %d \n',finclusion_total));
    
    %%%%%%%%%%
    %
    %   This section calculates the force on node 2
    %
    %%%%%%%%%%


    finclusion_2=zeros(1,3);
    
    finclusion_2(1) =  -(d/L_j)*finclusion_total(1) - const_general*((1/(2*L_j))*( ...
    -12*aI*B*b2*H_j^2*rhoI*z_1^2 - 16*B*b3*cI*H_j*rhoI*z_1^3 + 3*B*b2*eI*rhoI*z_1^4 - ...
	2*(10*G - 3*B*R_i*rhoI)*sqrt(H_j^2 + z_1^2)*(3*b3*cI*H_j*z_1 + 2*b2*(2*(-3 + dI)*H_j^2 - dI*z_1^2)) + ...
	12*aI*B*b2*H_j^2*rhoI*z_2^2 + 16*B*b3*cI*H_j*rhoI*z_2^3 - 3*B*b2*eI*rhoI*z_2^4 + ... 
	2*(10*G - 3*B*R_i*rhoI)*sqrt(H_j^2 + z_2^2)*(3*b3*cI*H_j*z_2 + 2*b2*(2*(-3 + dI)*H_j^2 - dI*z_2^2)) + ...
	6*b3*cI*H_j^3*(10*G - 3*B*R_i*rhoI)*log(z_1 + sqrt(H_j^2 + z_1^2)) - ...
	6*b3*cI*H_j^3*(10*G - 3*B*R_i*rhoI)*log(z_2 + sqrt(H_j^2 + z_2^2))));
    
    finclusion_2(2) = -(d/L_j)*finclusion_total(2) - const_general*(b1/(2*L_j)*(-4*dI*(10*G - 3*B*R_i*rhoI)*(H_j^2 + z_1^2)^(3/2) - ...
	3*B*eI*rhoI*(H_j^2 + z_1^2)^2 + 4*dI*(10*G - 3*B*R_i*rhoI)*(H_j^2 + z_2^2)^(3/2) + 3*B*eI*rhoI*(H_j^2 + z_2^2)^2));

    finclusion_2(3) = 0;

    %disp(sprintf('finclusion_2 outside=%d %d %d \n',finclusion_2));
    
    %%%%%%%%%%
    %
    %   This section calculates the force on node 1
    %
    %%%%%%%%%%

    finclusion_1=finclusion_total - finclusion_2;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [finclusion_1, finclusion_2]=ForceOutsideInclusion(MU,NU,B_2,R_i,MU_i,NU_i,B_i,b1_newbase,x1_newbase,x2_newbase,d,H_j,L_j,Rs,c,beta,YoungModulus_2)


    
    b1=b1_newbase(1);
    b2=b1_newbase(2);
    b3=b1_newbase(3);
    
    x=x1_newbase(1);
    z_1=x1_newbase(3)
    z_2=x2_newbase(3)
    MU=YoungModulus_2;
    

    %%%%%%%%%%
    %
    %   This section calculates the total force on the segment due to the inclusion
    %
    %%%%%%%%%%

    
    finclusion_total=zeros(1,3);
    finclusion_2=zeros(1,3);
    

    finclusion_total(1)=-MU*R_i^2*c*((3*H_j^4*b3 - 6*H_j^3*b2*z_1*sqrt((H_j^2 + z_1^2)/H_j^2) + 2*H_j^2*b2*sqrt(H_j^2 + z_1^2)*abs(z_1) + 3*H_j^2*b3*z_1^2 ...
        - 4*H_j*b2*z_1^3*sqrt((H_j^2 + z_1^2)/H_j^2) + 2*b2*z_1^2*sqrt(H_j^2 + z_1^2)*abs(z_1))*(H_j^4*NU + H_j^4 + 2*H_j^2*NU*z_2^2 + 2*H_j^2*z_2^2 + NU*z_2^4 + z_2^4)...
        - (3*H_j^4*b3 - 6*H_j^3*b2*z_2*sqrt((H_j^2 + z_2^2)/H_j^2) + 2*H_j^2*b2*sqrt(H_j^2 + z_2^2)*abs(z_2) + 3*H_j^2*b3*z_2^2 ...
        - 4*H_j*b2*z_2^3*sqrt((H_j^2 + z_2^2)/H_j^2) + 2*b2*z_2^2*sqrt(H_j^2 + z_2^2)*abs(z_2))*(H_j^4*NU + H_j^4 + 2*H_j^2*NU*z_1^2 ...
        + 2*H_j^2*z_1^2 + NU*z_1^4 + z_1^4))/(2*H_j^2*(H_j^4*NU + H_j^4 + 2*H_j^2*NU*z_1^2 + 2*H_j^2*z_1^2 + NU*z_1^4 + z_1^4)*(H_j^4*NU + H_j^4 + 2*H_j^2*NU*z_2^2 ...
        + 2*H_j^2*z_2^2 + NU*z_2^4 + z_2^4));

    finclusion_total(2)=MU*R_i^2*b1*c*(-sqrt(H_j^2 + z_1^2)*abs(z_2) + sqrt(H_j^2 + z_2^2)*abs(z_1))/(H_j^2*sqrt(H_j^2 + z_1^2)*sqrt(H_j^2 + z_2^2)*(NU + 1));

    finclusion_total(3)=0;

   % disp(sprintf('finclusion_total outside=%d %d %d \n',finclusion_total));

    %%%%%%%%%%
    %
    %   This section calculates the force on node 2
    %
    %%%%%%%%%%

    finclusion_2(1)=MU*R_i^2*c*(H_j^2*((H_j^2 + z_1^2)^2*(H_j^2*NU + H_j^2 + z_1^2*(NU + 1))*(4*H_j^2*b2*sqrt(H_j^2 + z_2^2)*(H_j^2*NU + H_j^2 + z_2^2*(NU + 1))*abs(H_j) ...
        - 4*b2*(H_j^2 + z_2^2)^(3/2)*(H_j^2*NU + H_j^2 + z_2^2*(NU + 1))*abs(H_j) + 3*b3*(H_j^2 + z_2^2)^2*(NU + 1)*(2*z_2*abs(H_j) ...
        + sqrt(-1/(NU + 1)^2)*(log((-H_j^2*NU*sqrt(-1/(NU + 1)^2) - H_j^2*sqrt(-1/(NU + 1)^2) + z_2*abs(H_j))/abs(H_j)) - log((H_j^2*NU*sqrt(-1/(NU + 1)^2) ...
        + H_j^2*sqrt(-1/(NU + 1)^2) + z_2*abs(H_j))/abs(H_j)))*(H_j^2*NU + H_j^2 + z_2^2*(NU + 1)))) - (H_j^2 + z_2^2)^2*(H_j^2*NU + H_j^2 ...
        + z_2^2*(NU + 1))*(4*H_j^2*b2*sqrt(H_j^2 + z_1^2)*(H_j^2*NU + H_j^2 + z_1^2*(NU + 1))*abs(H_j) - 4*b2*(H_j^2 + z_1^2)^(3/2)*(H_j^2*NU + H_j^2 + z_1^2*(NU + 1))*abs(H_j) ...
        + 3*b3*(H_j^2 + z_1^2)^2*(NU + 1)*(2*z_1*abs(H_j) + sqrt(-1/(NU + 1)^2)*(log((-H_j^2*NU*sqrt(-1/(NU + 1)^2) - H_j^2*sqrt(-1/(NU + 1)^2) + z_1*abs(H_j))/abs(H_j)) ...
        - log((H_j^2*NU*sqrt(-1/(NU + 1)^2) + H_j^2*sqrt(-1/(NU + 1)^2) + z_1*abs(H_j))/abs(H_j)))*(H_j^2*NU + H_j^2 + z_1^2*(NU + 1)))))*(H_j^4*NU + H_j^4 ...
        + 2*H_j^2*NU*z_1^2 + 2*H_j^2*z_1^2 + NU*z_1^4 + z_1^4)*(H_j^4*NU + H_j^4 + 2*H_j^2*NU*z_2^2 + 2*H_j^2*z_2^2 + NU*z_2^4 + z_2^4)*abs(H_j) ...
        + 2*H_j^2*d*(H_j^2 + z_1^2)^2*(H_j^2 + z_2^2)^2*(NU + 1)*((3*H_j^4*b3 - 6*H_j^3*b2*z_1*sqrt((H_j^2 + z_1^2)/H_j^2) + 2*H_j^2*b2*sqrt(H_j^2 + z_1^2)*abs(z_1) ...
        + 3*H_j^2*b3*z_1^2 - 4*H_j*b2*z_1^3*sqrt((H_j^2 + z_1^2)/H_j^2) + 2*b2*z_1^2*sqrt(H_j^2 + z_1^2)*abs(z_1))*(H_j^4*NU + H_j^4 + 2*H_j^2*NU*z_2^2 + 2*H_j^2*z_2^2 ...
        + NU*z_2^4 + z_2^4) - (3*H_j^4*b3 - 6*H_j^3*b2*z_2*sqrt((H_j^2 + z_2^2)/H_j^2) + 2*H_j^2*b2*sqrt(H_j^2 + z_2^2)*abs(z_2) + 3*H_j^2*b3*z_2^2 ...
        - 4*H_j*b2*z_2^3*sqrt((H_j^2 + z_2^2)/H_j^2) + 2*b2*z_2^2*sqrt(H_j^2 + z_2^2)*abs(z_2))*(H_j^4*NU + H_j^4 + 2*H_j^2*NU*z_1^2 + 2*H_j^2*z_1^2 ...
        + NU*z_1^4 + z_1^4))*(H_j^2*NU + H_j^2 + z_1^2*(NU + 1))*(H_j^2*NU + H_j^2 + z_2^2*(NU + 1)))/(4*H_j^2*H_j^2*L_j*(H_j^2 + z_1^2)^2*(H_j^2 + z_2^2)^2*(NU + 1)*(H_j^2*NU ...
        + H_j^2 + z_1^2*(NU + 1))*(H_j^2*NU + H_j^2 + z_2^2*(NU + 1))*(H_j^4*NU + H_j^4 + 2*H_j^2*NU*z_1^2 + 2*H_j^2*z_1^2 + NU*z_1^4 + z_1^4)*(H_j^4*NU ...
        + H_j^4 + 2*H_j^2*NU*z_2^2 + 2*H_j^2*z_2^2 + NU*z_2^4 + z_2^4));

    finclusion_2(2)=MU*R_i^2*b1*c*(H_j^2*(sqrt(H_j^2 + z_1^2) - sqrt(H_j^2 + z_2^2)) + d*(sqrt(H_j^2 + z_1^2)*abs(z_2) ...
        - sqrt(H_j^2 + z_2^2)*abs(z_1)))/(H_j^2*L_j*sqrt(H_j^2 + z_1^2)*sqrt(H_j^2 + z_2^2)*(NU + 1));

    finclusion_2(3)=0;

    disp(sprintf('finclusion_2 outside=%d %d %d \n',finclusion_2));
    
    %%%%%%%%%%
    %
    %   This section calculates the force on node 1
    %
    %%%%%%%%%%

    finclusion_1=finclusion_total - finclusion_2;
    disp(sprintf('finclusion_1 outside=%g %g %g \n',finclusion_1));
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segments=constructsegmentlist(rn,links)
[LINKMAX,m]=size(links);

segments=zeros(LINKMAX,14);
nseg=0;
for i=1:LINKMAX,
    n0=links(i,1);
    n1=links(i,2);
    if((n0~=0)&(n1~=0))
        nseg=nseg+1;
        segments(nseg,:)=[links(i,1:5),rn(n0,1:3),rn(n1,1:3),links(i,6:8)];
    end
end
segments=segments(1:nseg,:);

