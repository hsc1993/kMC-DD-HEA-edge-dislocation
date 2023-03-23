function links = generateLinks(rn,b,n)
    links = [];
    for i = 1:length(rn)-1
        links = [links;[i i+1 b n]];
    end
end
