%SCript to plot a pressure gradient plot and extract all of the required
%climb info

if ~exist('inclusion_pos_rad')
    inclusion_pos_rad=[0 10 0 5 0]; % 11 -11
end


x = (-15:0.01:15);
y = (-5:0.01:25);

for indx = 1:length(x)
    for indy = 1:length(y)
        
        if ((x(indx)-inclusion_pos_rad(1))^2 + (y(indy)-inclusion_pos_rad(2))^2 )^(0.5) < inclusion_pos_rad(4)
            s(indx,indy)=0;
        else
            s(indx,indy) = precipitatePressure([x(indx) y(indy) 0], inclusion_pos_rad);
        end
        
    end
end

imagesc(s)
set(gca, 'YDir', 'normal')
colormap default
h = colorbar;
h.TickLabels = 0:2:16;
set(h,'FontSize', 18)
t1 = ylabel(h, '\sigma_{rr} [a.u.]','FontSize', 25);
set(t1, 'FontSize', 22);


set(gca, 'xtick', []);
set(gca, 'ytick', []);
%xticklabels({'', '', '0', '5', '10', '15'})
hold on

%Normalize the emission data
[counts,centers] = hist(c1);


countn = counts /max(counts);
centern = centers/max(centers);


centerscale = centern*length(x)/3 + length(x)*6.3/10;
countscale = countn*length(x)/3 + length(x)/2;


plot([length(x)/2 max(centerscale)], [length(y)/2 , length(y)/2], 'r--', 'Linewidth', 2)
plot([length(x)/2 length(x)/2], [length(y)/2 , length(y)*7.3/8], 'r--', 'Linewidth', 2)



plot(centerscale, countscale,'k','Linewidth', 3)
% 
% xlabel('Distance From Precipitate Center [b]')
% ylabel('Vacancies Emitted [#]')

%Draw the x and y labels for the inlet graph
text(length(x)*2.05/3, length(y)*.95/2, 'Distance from Center', 'FontSize', 20)

y_label=text(length(x)*0.96/2, length(x)*2.01/3, 'Vacancies Emitted','FontSize', 20);
set(y_label, 'Rotation', 90);

%Now plot a sudo dislocation

%y1 = (0:250:3000);
% x1 = ones(length(y),1)*1200;

y1 = [0,250,700,900,1000,1250,1500,1700,2000,2250,2500,2750,3000];

x1=[1400;1400;1380;1180;1000;830;790;830;975;1280;1370;1400;1400];

scatter(x1,y1,100, [0.6350, 0.0780, 0.1840], 'filled');
plot(x1,y1 ,'color',[0.6350, 0.0780, 0.1840], 'Linewidth', 6);

%Label the dislocation
text(length(x)*2.38/4, length(y)*0.44/2, 'Dislocation Line', 'FontSize', 20);
annotation('textarrow', [0.46, 0.35], [0.23, .25]);


%Plot the glide plane
plot([length(x)/2 length(x)/2], [length(y)/2 , -length(y)*7.3/8], 'k--', 'Linewidth', 1)

%Label the glide plane
text(length(x)*2.46/4, length(y)*0.33/2, 'Glide Plane', 'FontSize', 20);
annotation('textarrow', [0.48, 0.38], [0.18, .20]);

%Adjust the axis
axis([length(x)/4.1 length(x)*.99 length(x)/11 length(x)*.93])