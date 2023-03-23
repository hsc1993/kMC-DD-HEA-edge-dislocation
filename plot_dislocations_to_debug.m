close all



check_id = 45;
hold on
plot3(rn_afterRelaxation(:,3),rn_afterRelaxation(:,1),rn_afterRelaxation(:,2),'x')
plot3(rn(:,3),rn(:,1),rn(:,2),'o')
scatter3(rn(check_id,3),rn(check_id,1),rn(check_id,2),'+')
hold off
legend('after relaxation','final rn')



zoom=5;
%view([-45,-45]); 
view([90,0]); 
% xlim([-plim/zoom*lattice_constant plim/zoom*lattice_constant]); ylim([-plim/zoom*0.5*lattice_constant plim/zoom*1.5*lattice_constant]); zlim([-plim/zoom*lattice_constant plim/zoom*lattice_constant]);
material shiny; %camlight;
%view(3); axis equal;
axis equal
%axis off
grid off

% view only the middle half, to exclude boundary condition
view(viewangle); 
portion = 1;
% xlim([-plim_x*portion plim_x*portion]); 
% ylim([0 plim_y*portion]); 
% zlim([-plim_z*portion plim_z*portion]);

grid on
daspect([1 1 1])




