function [stressMat] = globalStress(sigma,links,rn,b_vec, p)
%Function to calculate the global stress at a given point, p, given the
%sigma = applied uniform stress
%links = node connectivity
%rn = node positions


if isempty(sigma)
    stressMat = zeros(3,3); %Empty matrix if no applied stress
else
    stressMat = sigma; %Start off by adding the uniform applied stress
end


for i=1:length(links) %Loop through the links and add the stress
    
    %Get the node numbers
    node1 = links(i,1);
    node2 = links(i,2);
    
    %Get the node positions
    node1pos(1) = rn(node1,1);
    node1pos(2) = rn(node1,2);
    node1pos(3) = rn(node1,3);

    node2pos(1) = rn(node2,1);
    node2pos(2) = rn(node2,2);
    node2pos(3) = rn(node2,3);
    
    %Calcualte the local stress 
    temp_stress = stressStraight(node1pos,node2pos, b_vec, p);
    
    %Add it to the total stress
    stressMat=stressMat + temp_stress;
end


end

