function tau_rss_seglist = calRSS(fseg,rn,links,n,b_norm)

tau_rss_seglist = [];
for idx_seg = 1:size(links,1)
    idx_rn1 = links(idx_seg,1);
    idx_rn2 = links(idx_seg,2);
    l = norm(rn(idx_rn2,1:3)-rn(idx_rn1,1:3));
    t = (rn(idx_rn2,1:3)-rn(idx_rn1,1:3))/l; % line direction unit vector
    
    g = cross(n,t);  % glide direction unit vector
    tau_rss_seglist = [tau_rss_seglist;0.5*(fseg(idx_seg,1:3)+fseg(idx_seg,4:6))*g'/(b_norm*l)];
%       tau_rss_seglist in unites of [N/A^2]
end

end