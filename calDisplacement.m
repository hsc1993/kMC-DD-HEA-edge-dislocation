function totalDisplacement_current = calDisplacement(rn)
    temp = 0;

    for i = 1:size(rn,2)
        temp = temp+rn(i,1);
    end
    totalDisplacement_current = temp/size(rn,2);
end

