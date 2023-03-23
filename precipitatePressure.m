function [pressure] = precipitatePressure(pos,inclusion_pos_rad)
%Function to calculate the pressure felt at a position given information
%about the precipitate location and size

R = inclusion_pos_rad(4); %Radius of the inclusion

%Distance to the center of the precipitate
r = ( (pos(1)-inclusion_pos_rad(1))^2 + (pos(2)-inclusion_pos_rad(2))^2 + (pos(3)-inclusion_pos_rad(3))^2 )^0.5;

%Stand-in variables for E and m
E = 10E3; %Modulus
m = 8E4;

pressure = E*m/r; 

% if ~exist('NU')
%     NU=0.336;
% end
% 
% sigma_rr = -2*E*(R^2)* (m*R) /( r^3 * (NU+1));
% 
% sigma_theta = E*(R^2)* (m*R) /( r^3 * (NU+1));
% 
% %sigma_phi = E*(R^2)* (m*R) /( r^3 * (NU+1)); %Same thing as sigma theta
% 
% 
% %Reduced stress matrix to only consider the active portions of the pressure
% stress = [sigma_theta   0   0;
%           0 sigma_rr 0;
%           0 0   sigma_rr];
% 
% 
% %Calculate the pressure at the given point
% pressure = (stress(1,1) + stress(2,2) + stress(3,3))/3;
      
end

