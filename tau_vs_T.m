clear;close all;clc
Tlist = 100:100:1500;
tau_edge_list = [];
delta_tau_list = [];
Bedge_vgl_list = [];
Bedge_qm_list = [];

tau_total_list = [];
    





for i = 1:size(Tlist,2)
    T = Tlist(i)

    k=8.617e-5; % eV K-1
    alpha = 0.5;
    mu = 94; %GPa
    a0 = 3.24;%A
    b=sqrt(3)/2*a0;%A
    w=sqrt(6)/3*a0;%A
    
%     if T>=500
%         Ef = 0.13;
%     else
%         Ef = 0.02;
%     end

    Ef=0.7
%     Ef = 1.375

    rho = 5e12*1e-20;%A^-2
    mu0_prime = 1e13;%s^-1
    Em=0.32;%eV

%     Em=0.19;%eV
    
    dot_epsilon0 = 1e-04; %s^-1
    Bedge = 2.12e-10*1e-3;%GPa*s
    vgl = dot_epsilon0/rho/b;
    qm = mu0_prime*exp(-Em/k/T);
    
    tau_edge = (423.6-0.1*T)*1e-3; %GPa
    delta_tau = alpha*mu*b/w*exp(-Ef/k/T);%GPa
    tau_glide = Bedge*vgl/b;
    tau_jump = -2*Bedge*qm;
    tau_total= tau_edge+delta_tau+tau_glide+tau_jump;
    tau_edge_list = [tau_edge_list,tau_edge];
    delta_tau_list = [delta_tau_list,delta_tau];
    Bedge_vgl_list = [Bedge_vgl_list,tau_glide];
    Bedge_qm_list = [Bedge_qm_list,tau_jump];
    tau_total_list = [tau_total_list,tau_total];

end

subplot(2,3,1)
plot(Tlist,tau_edge_list,'LineWidth',4)
title('tau_{edge}', 'FontSize', 30)
xlabel('T [K]', 'FontSize', 30)
ylabel('stress [GPa]', 'FontSize', 30)

subplot(2,3,2)
plot(Tlist,delta_tau_list,'LineWidth',4)
title('\delta tau^*', 'FontSize', 30)
xlabel('T [K]', 'FontSize', 30)
ylabel('stress [GPa]', 'FontSize', 30)

subplot(2,3,3)
plot(Tlist,Bedge_vgl_list,'LineWidth',4)
title('tau_{glide}', 'FontSize', 30)
xlabel('T [K]', 'FontSize', 30)
ylabel('stress [GPa]', 'FontSize', 30)

subplot(2,3,4)
plot(Tlist,Bedge_qm_list,'LineWidth',4)
title('tau_{migration}', 'FontSize', 30)
xlabel('T [K]', 'FontSize', 30)
ylabel('stress [GPa]', 'FontSize', 30)

subplot(2,3,5)
plot(Tlist,tau_total_list,'LineWidth',4)
title('tau_{total}', 'FontSize', 30)
xlabel('T [K]', 'FontSize', 30)
ylabel('stress [GPa]', 'FontSize', 30)


