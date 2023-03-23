function typelist = generateNewTypelist(rn,links)
typelist = [];
    for idx_seg = 1:size(links,1)
        idx_r1 = links(idx_seg,1);
        idx_r2 = links(idx_seg,2);
        if (rn(idx_r1,4) == 0) && (rn(idx_r2,4)==0)
    %         straight
            typelist = [typelist;1];
        end
        if (rn(idx_r1,4) == 7) && (rn(idx_r2,4)==7) && abs(rn(idx_r1,2)-rn(idx_r2,2))<1
    %         small straight line between two jogs
            typelist = [typelist;4];
        end
        if (rn(idx_r1,4) == 7) && (rn(idx_r2,4)==7) && abs(rn(idx_r1,2)-rn(idx_r2,2))>1
    %         jog
            typelist = [typelist;2];
        end

        if (rn(idx_r1,4) == 0) && (rn(idx_r2,4)==7)
    %         connection
            typelist = [typelist;3];
        end
        if (rn(idx_r1,4) == 7) && (rn(idx_r2,4)==0)
    %         connection
            typelist = [typelist;3];
        end
    end
end

