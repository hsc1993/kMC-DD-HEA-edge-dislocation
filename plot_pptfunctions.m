clear all

steps = 154;
step = 0.2;
finclusion_1=zeros(steps,1);
finclusion_2=zeros(steps,1);
fit1=zeros(steps,1);
fit2=zeros(steps,1);
fi21=zeros(steps,1);
fi22=zeros(steps,1);
f=zeros(steps,1);
EO = 1;
beta = 28.35;
c = 0.17;
R =2.33;
vO = 0.3;
Rs = 33;
b1=1;
b2=-1;
b3=0;
h1=-0.5;
h2=0.5;


for i = 1:steps
    radius(i)=R+(i-1)*step;
    r = radius(i);
    H=sqrt(r^2-h1^2);
    z=h1; %sqrt(r^2-H^2);
    z_1=h1;
    z_2=h2;
    d=h2;
    L = h2-h1;
    
    u(i)=c*(1-(r-R)/(Rs-R)*exp((Rs-r)/beta));
    
    err(i) = (c*exp((-r + Rs)/beta)*(beta - r + R))/(beta*(-R + Rs));
    
    ett(i) = u(i)/r;
    
    sorr1(i)=EO/((1+vO)*(1-2*vO))*((1-vO)*err(i)+2*vO*ett(i));
    sott1(i)=EO/((1+vO)*(1-2*vO))*(ett(i)+vO*err(i));
    
    sorr2(i) = EO*c*(2*vO*beta*(-R + Rs + (R - r)*exp((Rs - r)/beta)) + r*(vO - 1)*(R + beta - r)*exp((Rs - r)/beta))/(beta*r*(vO + 1)*(2*vO - 1)*(R - Rs));
    sott2(i) = -EO*c*(vO*r*(R + beta - r)*exp((Rs - r)/beta) - beta*(-R + Rs + (R - r)*exp((Rs - r)/beta)))/(beta*r*(vO + 1)*(2*vO - 1)*(R - Rs));
    
    sorr3(i)=-((c*exp((-r + Rs)/beta)*EO*(r*(r - R)*(-1 + vO) + beta*(r + r*vO - 2*R*vO)))/ ...
    (beta*r*(-R + Rs)*(1 + vO)*(-1 + 2*vO)));

    sott3(i)=-((c*exp((-r + Rs)/beta)*EO*(r*(-r + R)*vO + beta*(r - R + r*vO)))/ ...
    (beta*r*(-R + Rs)*(1 + vO)*(-1 + 2*vO)));

	ft11 = (b2*c*exp((Rs - sqrt(H^2 + z^2))/beta)*EO*H^2*(beta*R*(beta + sqrt(H^2 + z^2)) + ...
       exp(sqrt(H^2 + z^2)/beta)*(2*beta + R)*(H^2 + z^2)*real(-expint(-(sqrt(H^2 + z^2)/ ...
          beta))))/(2*beta^2*(-R + Rs)*(1 + vO)*(H^2 + z^2)));
	  
    ft12 = (b2*c*exp((Rs -sqrt(H^2 + z^2))/beta)*EO*(R + beta*vO - sqrt(H^2 + z^2) + ... 
       exp(sqrt(H^2 + z^2)/beta)* ...
         R*vO*real(-expint(-(sqrt(H^2 + z^2)/beta))))/((-R + Rs)* ...
     (1 + vO)*(-1 + 2*vO)));

    ft13 = -((b3*c*exp((Rs - sqrt(H^2 + z^2))/beta)*EO*H* ...
    sqrt(z^2)*(-R + sqrt(H^2 + z^2)))/((-R + Rs)*(1 + vO)*z*sqrt(H^2 + z^2)));
  
    ft2 = -((b1*c*exp((Rs - sqrt(H^2 + z^2))/beta)* ...
       EO*(R + beta*vO - sqrt(H^2 + z^2) + ...
        exp(sqrt(H^2 + z^2)/beta)* ...
          R*vO*real(-expint(-(sqrt(H^2 + z^2)/beta))))/((-R + Rs)*(1 + ...
        vO)*(-1 + 2*vO))));
    f1(i)=ft11+ft12+ft13;
    f2(i)=ft2;
    
    const_general=c*exp(Rs/beta)*EO/((-R + Rs)*(1 + vO));
    
    f_total_1(1)= b2*H^2/(2*beta^2)*(-((exp(-(sqrt (H^2 + z_1^2)/beta))*(beta*R*(beta + sqrt(H^2 + z_1^2)) + ...
    exp(sqrt (H^2 + z_1^2)/beta)*(2*beta + R)*(H^2 + z_1^2)*real(-expint(-(-(sqrt (H^2 + z_1^2)/beta))))))/(H^2 + z_1^2)) + ...
    (exp(-(sqrt (H^2 + z_2^2)/beta))*(beta*R*(beta + sqrt(H^2 + z_2^2)) + ...
    exp(sqrt (H^2 + z_2^2)/beta)*(2*beta + R)*(H^2 + z_2^2)*real(-expint(-(-(sqrt (H^2 + z_2^2)/beta))))))/(H^2 + z_2^2));

    f_total_2(1)=b1/(-1 + 2*vO)*(exp(-(sqrt (H^2 + z_1^2)/beta))*(R + beta*vO - sqrt(H^2 + z_1^2) + ...
    exp(sqrt (H^2 + z_1^2)/beta)*R*vO*real(-expint(-(-(sqrt (H^2 + z_1^2)/beta))))) - ...
    exp(-(sqrt (H^2 + z_2^2)/beta))*(R + beta*vO - sqrt(H^2 + z_2^2) + ...
    exp(sqrt (H^2 + z_2^2)/beta)*R*vO*real(-expint(-(-(sqrt (H^2 + z_2^2)/beta))))));

    f_total_3(1)=b3*H*((exp(-(sqrt(H^2 + z_1^2)/beta))*sqrt(z_1^2)*(R - sqrt(H^2 + z_1^2)))/(z_1*sqrt(H^2 + z_1^2)) + ...
    (exp(-(sqrt (H^2 + z_2^2)/beta))*sqrt(z_2^2)*(-R + sqrt(H^2 + z_2^2)))/(z_2*sqrt(H^2 + z_2^2)));

    fit1(i)=const_general*(f_total_1(1) + f_total_2(1) + f_total_3(1));

    fit2(i)=-const_general*(b1/(-1 + 2*vO)*(exp(-(sqrt(H^2 + z_1^2)/beta))*(R + beta*vO - sqrt(H^2 + z_1^2) + ...
    (exp((sqrt(H^2 + z_1^2)/beta))*R*vO*real(-expint(-(-(sqrt(H^2 + z_1^2)/beta))))) - ...
    (exp(-(sqrt(H^2 + z_2^2)/beta))*(R + beta*vO - sqrt(H^2 + z_2^2) + ...
    (exp((sqrt(H^2 + z_2^2)/beta))*R*vO*real(-expint(-(-(sqrt(H^2 + z_2^2)/beta))))))))));

    f_2_1(1) = -(d/L)*f_total_1(1) - ((b2*H^2)/(L)*(exp(-(sqrt(H^2 + z_1^2)/beta))*(-1 + R/sqrt(H^2 + z_1^2)) + ...
    exp(-(sqrt(H^2 + z_2^2)/beta))*(1 - R/sqrt(H^2 + z_2^2))));

    f_2_2(1) = -(d/L)*f_total_2(1) - (b2/(L*(-1 + 2*vO))*(exp(-(sqrt(H^2 + z_1^2)/beta))*(H^2 - beta*(beta - R)*(-1 + vO) + ...
    z_1^2 - (R + beta*(-1 + vO))*sqrt(H^2 + z_1^2)) + ...
    exp(-(sqrt(H^2 + z_2^2)/beta))*(-H^2 + beta*(beta - R)*(-1 + vO) - ...
    z_2^2 + (R + beta*(-1 + vO))*sqrt(H^2 + z_2^2))));

    f_2_3(1)=-(d/L)*f_total_3(1) - (b3*H/(2*L*beta^2)*(-( ...
    (exp(-(sqrt(H^2 + z_1^2)/beta))*(-beta*(R*sqrt(H^2 + z_1^2)*H^2 +...
    2*beta^2*(H^2 + z_1^2) + beta*(-H^2*R - 2*R*z_1^2 + 2*(H^2 + z_1^2)^(3/2))) - ...
    exp(sqrt(H^2 + z_1^2)/beta)*(2*beta*H^2 + 2*beta^2*R + H^2*R)*(H^2 + z_1^2) ...
    *real(-expint(-(-(sqrt(H^2 + z_1^2)/beta))))))/(H^2 + z_1^2)) + ...
    (exp(-(sqrt(H^2 + z_2^2)/beta))*(-beta*(H^2*R*sqrt(H^2 + z_2^2) + ...
    2*beta^2*(H^2 + z_2^2) + beta*(-H^2*R - 2*R*z_2^2 + 2*(H^2 + z_2^2)^(3/2))) - ... 
    exp(sqrt(H^2 + z_2^2)/beta)*(2*beta*H^2 + 2*beta^2*R + H^2*R)*(H^2 + z_2^2) ...
    *real(-expint(-(-(sqrt(H^2 + z_2^2)/beta))))))/(H^2 + z_2^2)));

    fi21(i)=const_general*(f_2_1(1) + f_2_2(1) + f_2_3(1));

    fi22(i)=-(d/L)*ft2 + const_general*(b1/(L*(-1 + 2*vO))*(exp(-(sqrt(H^2 + z_1^2)/beta))*(-H^2 + beta*(beta - R)*(-1 + vO) - ...
    z_1^2 + (R + beta*(-1 + vO))*sqrt(H^2 + z_1^2)) + exp(-(sqrt(H^2 + z_2^2)/beta))*(H^2 - beta*(beta - R)*(-1 + vO) + ... 
    z_2^2 - (R + beta*(-1 + vO))*sqrt(H^2 + z_2^2))));

    finclusion_1(i)=fit1(i) - fi21(i);
    finclusion_2(i)=fit2(i) - fi22(i);
    f(i)=norm([finclusion_1(i),finclusion_2(i)]);
    
end

%plot(radius,fi21);
%hold on
%plot(radius,fi22);
%plot(radius,fit1);
%plot(radius,fit2);
%plot(radius,finclusion_1);
%plot(radius,finclusion_2);
%plot(radius,f);
%hold off
%legend('fi21','fi22','fit1','fit2','finclusion 1','finclusion 2','f');

plot(radius,sorr1);
hold on
plot(radius,sott1);
plot(radius,sorr2);
plot(radius,sott2);
plot(radius,sorr3);
plot(radius,sott3);
hold off
legend('\sigma_{rr} 1', '\sigma_{\theta\theta} 1','\sigma_{rr} 2', '\sigma_{\theta\theta} 2','\sigma_{rr} 3', '\sigma_{\theta\theta} 3');
xlabel('Radius r');
ylabel('Stress');

% plot(radius,err);
% hold on
% plot(radius,ett);
% hold off
% legend('\epsilon_{rr}', '\epsilon_{\theta\theta}');
% xlabel('Radius r');
% ylabel('Strain');

%plot(radius,u);
%xlabel('Radius r');
%ylabel('Displacement');
