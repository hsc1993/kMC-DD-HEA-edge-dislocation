function [stress] = precipitateStress(pos,inclusion_pos_rad)
%Function to calculate the stress felt at a position given information
%about the precipitate location and size

R = inclusion_pos_rad(4); %Radius of the inclusion

%Distance to the center of the precipitate
r = ( (pos(1)-inclusion_pos_rad(1))^2 + (pos(2)-inclusion_pos_rad(2))^2 + (pos(3)-inclusion_pos_rad(3))^2 )^0.5;

%Stand-in variables for E and m
E = 10E3; %Modulus
m = 0.2;

if ~exist('NU')
    NU=0.336;
end

sigma_rr = -2*E*(R^2)* (m*R) /( r^3 * (NU+1));

sigma_theta = E*(R^2)* (m*R) /( r^3 * (NU+1));

%sigma_phi = E*(R^2)* (m*R) /( r^3 * (NU+1)); %Same thing as sigma theta

%Construct a vector from the center of the inclusion to the point of
%interest
vec(1) = pos(1) - inclusion_pos_rad(1);
vec(2) = pos(2) - inclusion_pos_rad(2);
vec(3) = pos(3) - inclusion_pos_rad(3);

vectical_vec = [0 0 1];

%Calculate the vectical angle to use the simplified stress format
theta = acosd(dot(vec,vectical_vec)./ (norm(vec)*norm(vectical_vec)));

stress = [sigma_theta   0   0;
          0 (sigma_rr*sind(theta)^2 + sigma_theta*cosd(theta)^2) (sigma_rr*sind(theta)*cosd(theta) + sigma_theta*cosd(theta)*sind(theta));
          0 (sigma_rr*sind(theta)*cosd(theta) + sigma_theta*cosd(theta)*sind(theta))    (sigma_rr*sind(theta)^2 + sigma_theta*cosd(theta)^2)];


%Check to see if there is any rotation about the z-axis
y_vec = [1 0 0];

%Calculate the vectical angle to use the simplified stress format
phi = acosd(dot(vec,y_vec)./ (norm(vec)*norm(y_vec)));

%Only perform transformation if the angle is not 90
if (phi~=90)
    
   %Transformation matrix
   B = [0   0   1;
       cosd(phi)    sind(phi)  0;
       -sind(phi)   cosd(phi)   0];
   
   
   stress = transpose(B)*stress*B;
   
end
      
%%
% 
% %Attempt #2 to deal with rotation matrix

% %Construct a vector from the center of the inclusion to the point of
% %interest
% new_x(1) = pos(1) - inclusion_pos_rad(1);
% new_x(2) = pos(2) - inclusion_pos_rad(2);
% new_x(3) = pos(3) - inclusion_pos_rad(3);
% 
% old_x = [1 0 0];
% old_y = [0 1 0];
% old_z = [0 0 1];
% 
% %Find the new y and z vectors by calculating the cross products
% new_y = cross(old_z, new_x);
% new_z = cross(new_x, new_y);
% 
% %Calculate the vectical angle to use the simplified stress format
% thetax = acosd(dot(new_x,old_x)./ (norm(new_x)*norm(old_x)));
% thetay = acosd(dot(new_y,old_y)./ (norm(new_y)*norm(old_y)));
% thetaz = acosd(dot(new_z,old_z)./ (norm(new_z)*norm(old_z)));
% 
% 
% Rx = [  1   0   0;
%         0 cosd(thetax)  -sind(thetax);
%         0   sind(thetax)    cosd(thetax)];
%     
% Ry = [cosd(thetay)  0   sind(thetay);
%       0 1   0;
%       -sind(thetay)         0   cosd(thetay)];
%   
% Rz = [cosd(thetaz)  -sind(thetaz)   0;
%       sind(thetaz)  cosd(thetaz)    0;
%       0 0   1];
% 
% 
% %Calculate the rotation matrix
% B = Rz*Ry*Rx;
% 
% stress = [sigma_rr  0   0;
%           0     sigma_theta     0;
%           0     0               sigma_theta];
%       
% %Rotate the stress to the desired x-y-z plane
% stress = transpose(B)*stress*B;

end

