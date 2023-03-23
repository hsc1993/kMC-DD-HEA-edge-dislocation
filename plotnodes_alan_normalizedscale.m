%function to plot dislocation segments
function plotnodes_alan(rn,links,plim,doinclusion,inclusion_pos_rad,SurfacePlane,typelist,viewangle,plim_x,plim_y,plim_z,includeSuperjogs)
%plot nodes
%only those nodes within [-plim,plim] in x,y,z directions are plotted
eps=1e-12;
perfectburg=1/2.*[1 1 0];
normshockley=norm(1/6.*[1 1 2]);
normstairrod=norm(1/6.*[1 1 0]);
normfrank=norm(1/3.*[1 1 1]);
normhirth=norm(1/3.*[1 0 0]);
plot3(0,0,0); hold on;
LINKMAX=length(links(:,1));
lattice_constant=0.4032;
lattice_constant=1;


for i=1:LINKMAX,
    n0=links(i,1);
    n1=links(i,2);
    b=links(i,3:5);
    if includeSuperjogs ==1 
%     if((n0~=0)&(n1~=0)&(max(max(abs(rn([n0,n1],:))))<=plim))
%             if (rn(n0,1) == rn(n1,1)) && (rn(n0,2) == rn(n1,2)) &&(rn(n0,3) == rn(n1,3))
%                 continue % two nodes being at the same position 
%             end

            if n0==size(rn,1) || n1==size(rn,1)
                continue % num seg = num node - 1, cant plot seg(maxIdxNode)
            end
            
            if typelist(n0) == 1 % straight seg
%                 h=plot3t(rn([n0,n1],1)*lattice_constant,rn([n0,n1],2)*lattice_constant,rn([n0,n1],3)*lattice_constant,lattice_constant/sqrt(2)/5,'b');
                h=plot3t(rn([n0,n1],3)*lattice_constant,rn([n0,n1],1)*lattice_constant,rn([n0,n1],2)*lattice_constant,lattice_constant/2,'b');
            elseif   typelist(n0) == 2 % jog
%                 h=plot3t(rn([n0,n1],1)*lattice_constant,rn([n0,n1],2)*lattice_constant,rn([n0,n1],3)*lattice_constant,lattice_constant/sqrt(2)/5,'r');
                h=plot3t(rn([n0,n1],3)*lattice_constant,rn([n0,n1],1)*lattice_constant,rn([n0,n1],2)*lattice_constant,lattice_constant/2,'b');
            else  % connection seg
%                 h=plot3t(rn([n0,n1],1)*lattice_constant,rn([n0,n1],2)*lattice_constant,rn([n0,n1],3)*lattice_constant,lattice_constant/sqrt(2)/5,'k');
                h=plot3t(rn([n0,n1],3)*lattice_constant,rn([n0,n1],1)*lattice_constant,rn([n0,n1],2)*lattice_constant,lattice_constant/2,'b');
            end
            hold on
%                 scatter3(rn([n0,n1],3)*lattice_constant,rn([n0,n1],1)*lattice_constant,rn([n0,n1],2)*lattice_constant);

%             plot direction for x y z
%             plot3([10 15],[0 0],[0 0],'LineWidth',4,'Color','r')
%             plot3([10 10],[0 5],[0 0],'LineWidth',4,'Color','b')
%             plot3([10 10],[0 0],[0 5],'LineWidth',4,'Color','k')
            hold off
            set(h, 'EdgeLighting','gouraud','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
            %material shiny;
    else
            h=plot3t(rn([n0,n1],3)*lattice_constant,rn([n0,n1],1)*lattice_constant,rn([n0,n1],2)*lattice_constant,lattice_constant/sqrt(2)/5,'b');

            set(h, 'EdgeLighting','gouraud','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
            %material shiny;
    end
end





% for i=1:LINKMAX,
%     n0=links(i,1);
%     n1=links(i,2);
%     b=links(i,3:5);
%     if((n0~=0)&(n1~=0)&(max(max(abs(rn([n0,n1],:))))<=plim))
%         %filter out "infinity" lines
%         if(abs(norm(b)-normshockley)<eps)
%             plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),'b.-');
%         elseif(abs((norm(b)-norm(perfectburg)))<eps)
%             plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),'r.-');
%         elseif(abs(norm(b)-normstairrod)<eps)
%             plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),'g.-');
%             %disp('A stair-rod has formed');
%         elseif(abs(norm(b)-normfrank)<eps)
%             plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),'k.-');
%         elseif(abs(norm(b)-normhirth)<eps)
%             plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),'y.-');
%         else
%             plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),'m.-');
%         end
%     end
% end
% xlabel('x - $$[111]$$ [nm]','interpreter','latex','FontSize',18,'FontWeight','bold');
% ylabel('y - $$[1\bar{1}0]$$ [nm]','interpreter','latex','FontSize',18,'FontWeight','bold');
% zlabel('z - $$[11\bar{2}]$$ [nm]','interpreter','latex','FontSize',18,'FontWeight','bold');
ylabel('x - $$[111]$$ [A]','interpreter','latex','FontSize',18,'FontWeight','bold');
zlabel('y - $$[1\bar{1}0]$$ [A]','interpreter','latex','FontSize',18,'FontWeight','bold');
xlabel('z - $$[11\bar{2}]$$ [A]','interpreter','latex','FontSize',18,'FontWeight','bold');
set(gca, 'FontSize', 16)
if(doinclusion~=0)
    for(i=1:size(inclusion_pos_rad,1))
        [x,y,z]=ellipsoid(inclusion_pos_rad(i,1)*lattice_constant,inclusion_pos_rad(i,2)*lattice_constant,inclusion_pos_rad(i,3)*lattice_constant,inclusion_pos_rad(i,4)*lattice_constant,inclusion_pos_rad(i,4)*lattice_constant,inclusion_pos_rad(i,4)*lattice_constant,12);
        surf(x,y,z,'EdgeColor','none','FaceColor','interp','FaceLighting','gouraud');
        alpha(0.9);
        colormap bone;
    end
end


%mark precipitate nodes
% for i=1:length(rn(:,4))
%     if rn(i,4) == 3
%         plot3(rn(i,1),rn(i,2),rn(i,3),'co');
%     end
%     if rn(i,4) == 2
%         plot3(rn(i,1),rn(i,2),rn(i,3),'go');
%     end
% end
% if(norm(SurfacePlane(1:3)) ~= 0)
%     if SurfacePlane(2) ~= 0
%         [X,Z] = meshgrid(-plim:500:plim);
%         Y = -(SurfacePlane(1) * X + SurfacePlane(3) * Z + SurfacePlane(4))/SurfacePlane(2);
%         h=surf(X,Y,Z);
%         set(h,'FaceColor',[1 0 0],'FaceAlpha',0.5,'EdgeAlpha',0);
%     elseif SurfacePlane(1) ~= 0
%         [Y,Z] = meshgrid(-plim:500:plim);
%         X = -(SurfacePlane(2) * Y + SurfacePlane(3) * Z + SurfacePlane(4))/SurfacePlane(1);
%         h=surf(X,Y,Z);
%         set(h,'FaceColor',[1 0 0],'FaceAlpha',0.5,'EdgeAlpha',0);
%     elseif SurfacePlane(3) ~= 0
%         [X,Y] = meshgrid(-plim:500:plim);
%         Z = -(SurfacePlane(1) * X + SurfacePlane(2) * Y + SurfacePlane(4))/SurfacePlane(3);
%         h=surf(X,Y,Z);
%         set(h,'FaceColor',[1 0 0],'FaceAlpha',0.5,'EdgeAlpha',0);
%     end
%     [Y,Z] = meshgrid(-plim:500:plim);
%     X = 0 * Y + 0 * Z - 2000;
%     h=surf(X,Y,Z);
%     set(h,'FaceColor',[0 1 0],'FaceAlpha',0.5,'EdgeAlpha',0);
% end
hold off
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
xlim([plim_x(1)*portion plim_x(2)*portion]); 
ylim([plim_y(1)*portion plim_y(2)*portion]); 
% ylim([-plim_y*portion plim_y*portion]); 
zlim([0 plim_z*portion]);

grid on
daspect([10 10 10])


