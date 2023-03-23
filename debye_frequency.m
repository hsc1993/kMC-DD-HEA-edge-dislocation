function [a_alloy,debye_frequency_alloy] = debye_frequency()



% physical constants
kb = 8.617333262145e-5;	% eV * K^-1
h_bar = 6.582119569e-16; % eV * s

% lattice constant in A
a_Nb = 3.30;
a_Mo = 3.15;
a_Ta = 3.31;
a_W = 3.16;
a_alloy = 0.25*(a_Nb+a_Mo+a_Ta+a_W);
edgeburger = a_alloy/2*sqrt(3);% A
screwburger = a_alloy/2*sqrt(3);% A

%% debye_frequency_alloy
debye_T_Nb = 275;
debye_T_Mo = 450;
debye_T_Ta = 240;
debye_T_W = 400;

debye_frequency_Nb = kb*debye_T_Nb/h_bar;  % s^-1
debye_frequency_Mo = kb*debye_T_Mo/h_bar;
debye_frequency_Ta = kb*debye_T_Ta/h_bar;
debye_frequency_W = kb*debye_T_W/h_bar;
debye_frequency_alloy = 0.25*(debye_frequency_Nb+debye_frequency_Mo+debye_frequency_Ta+debye_frequency_W);

