function [numbins,weights,energies,cumulative_weights] = Emig_data_histing()

fileID = fopen('Emig_ed.txt','r');
formatSpec = '%f';
sizeA = [1 Inf];
A = fscanf(fileID,formatSpec,sizeA);
edges = 0:0.1:0.8;

% h = histogram(A,edges);
% xticks(edges)

% numbins = h.NumBins
% weights = h.Values
% energies = [];
% cumulative_weights = [weights(1)]

numbins = 8;
weights = [12    18    26    13    19    11     4     3];
energies = [];
cumulative_weights = [12];

h_BinEdges = [0    0.1000    0.2000    0.3000    0.4000    0.5000    0.6000    0.7000    0.8000];

h_Values = [12    18    26    13    19    11     4     3];

for i = 1:numbins
    energies = [energies,0.5*(h_BinEdges(i)+h_BinEdges(i+1))];
end
for i = 2:numbins
    cumulative_weights = [cumulative_weights,cumulative_weights(i-1)+h_Values(i)];
end


end
