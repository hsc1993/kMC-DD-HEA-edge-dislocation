function [Evac,Evac_mean] = Evac_sampling(numbins,weights,energies,cumulative_weights)


samplingnum = randi([1,sum(weights)]);

sample_i = 1;

if samplingnum == cumulative_weights(end)
    sample_i = numbins;
else
    while cumulative_weights(sample_i)<=samplingnum
        sample_i = sample_i+1;
    end
end

Evac = energies(sample_i);


%% calculate Efv_mean
Evac_mean = 0;
for i = 1:numbins
    Evac_mean = Evac_mean+energies(i)*weights(i);
end
Evac_mean = Evac_mean/cumulative_weights(end);


end




