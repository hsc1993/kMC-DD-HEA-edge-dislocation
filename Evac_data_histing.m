function [numbins,weights,energies,cumulative_weights] = Evac_data_histing()

fileID = fopen('Evac_ed.txt','r');
formatSpec = '%f';
sizeA = [1 Inf];
A = fscanf(fileID,formatSpec,sizeA);
edges = -0.25:0.25:1.75;
% ax.subplot(2,2,1)

% h = histogram(A,edges);
% xticks(edges)

% numbins = h.NumBins
% weights = h.Values
% energies = [];
% cumulative_weights = [weights(1)]

numbins = 8;
weights = [1     1    16    17    42    37    24     1];
energies = [];
cumulative_weights = [1];

h_BinEdges = [-0.2500         0    0.2500    0.5000    0.7500    1.0000    1.2500    1.5000    1.7500];
h_Values = [1     1    16    17    42    37    24     1];
for i = 1:numbins
    energies = [energies,0.5*(h_BinEdges(i)+h_BinEdges(i+1))];
end
for i = 2:numbins
    cumulative_weights = [cumulative_weights,cumulative_weights(i-1)+h_Values(i)];
end


end
