function [ep_inc,wp_inc,L]=calcplasticstrainincrement(rnnew,rn,links_new,links,Volume)

len_new=size(links_new,1);
len=size(links,1);
%disp(sprintf('len new %d len %d \n',size(links,1),size(links_new,1)))

if len==len_new

    seg=rn(links(:,2),1:3)-rn(links(:,1),1:3);
    segnew=rnnew(links(:,2),1:3)- rnnew(links(:,1),1:3);
    
    dx1=rnnew(links(:,2),1:3)-rn(links(:,1),1:3); % Diagonal
    
    dA=cross(segnew+seg,dx1);

    fp_inc=0.5.*(links(:,3:5)'*dA)./Volume;
    ep_inc=0.5.*(fp_inc+fp_inc');
    wp_inc=0.5.*(fp_inc-fp_inc');
    
elseif len<len_new
    
    %disp(sprintf('Strain inc len < len_new %d %d\n',len,len_new));
    
    seg1=rn(links(:,2),1:3);
    seg2=rn(links(:,1),1:3);
    seg=seg1-seg2;

    
    segnew=rnnew(links(:,2),1:3)- rnnew(links(:,1),1:3);
    
    dx1=rnnew(links(:,2),1:3)-rn(links(:,1),1:3); % Diagonal
    
    dA=cross(segnew+seg,dx1);

    fp_inc=0.5.*(links(:,3:5)'*dA)./Volume;
    ep_inc=0.5.*(fp_inc+fp_inc');
    wp_inc=0.5.*(fp_inc-fp_inc');
    
elseif len>len_new
    
    %disp(sprintf('Strain inc len > len_new %d %d\n',len,len_new));
    
    seg=rn(links_new(:,2),1:3)-rn(links_new(:,1),1:3);
    segnew=rnnew(links_new(:,2),1:3)- rnnew(links_new(:,1),1:3);
    
    dx1=rnnew(links_new(:,2),1:3)-rn(links_new(:,1),1:3); % Diagonal
    
    dA=cross(segnew+seg,dx1);

    fp_inc=0.5.*(links_new(:,3:5)'*dA)./Volume;
    ep_inc=0.5.*(fp_inc+fp_inc');
    wp_inc=0.5.*(fp_inc-fp_inc');

end

L=sum(sqrt(sum((rnnew(links_new(:,2),1:3)- rnnew(links_new(:,1),1:3)).^2,2)),1);

% glideplane=[ -1 -1 -1 ]/norm([ -1 -1 -1 ]);
% tot_epsilon=(ep_inc*glideplane')';
% shear_strain=cross(cross(glideplane,tot_epsilon),glideplane);
% norm_shear_strain=norm(shear_strain)

%fp_inc=0.5.*(links(:,3:5)'*dA)./Volume;
%ep_inc=0.5.*(fp_inc+fp_inc');
%wp_inc=0.5.*(fp_inc-fp_inc');