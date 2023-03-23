function rn_afterDD = sortrn(rn_afterDD)
    rn_afterDD_new = rn_afterDD;

    rn_temp = sortrows(rn_afterDD,3);
    rn_temp_idx = find(rn_temp(:,4)==7);
    rn_temp_idx = rn_temp_idx(2:end-1);
    rn_temp = rn_temp(rn_temp_idx,:);

    numsuperjog = size(rn_temp,1)/4;
    for iter = 1:numsuperjog
        temp = 4*(iter-1)+1;
        idx_rn_superjog = [temp temp+1 temp+2 temp+3];
        rn_superjog = rn_temp(idx_rn_superjog,:);
        % check if 1,2 nodes need to be swapped
        if rn_superjog(1,2)>rn_superjog(2,2)
            swap_idx1 = rn_temp_idx(idx_rn_superjog(1));
            swap_idx2 = rn_temp_idx(idx_rn_superjog(2));

            swap_temp = rn_afterDD(swap_idx1,:);
            rn_afterDD_new(swap_idx1,:) = rn_afterDD_new(swap_idx2,:) ;
            rn_afterDD_new(swap_idx2,:) = swap_temp;
        end
        % check if 3,4 nodes need to be swapped
        if rn_superjog(3,2)<rn_superjog(4,2)
            swap_idx3 = rn_temp_idx(idx_rn_superjog(3));
            swap_idx4 = rn_temp_idx(idx_rn_superjog(4));

            swap_temp = rn_afterDD(swap_idx3,:);
            rn_afterDD_new(swap_idx3,:) = rn_afterDD_new(swap_idx4,:) ;
            rn_afterDD_new(swap_idx4,:) = swap_temp;
        end
    end

end