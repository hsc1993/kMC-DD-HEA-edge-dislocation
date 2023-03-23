function [Emig,Emig_mean] = Emig_sampling(numbins,weights,energies,cumulative_weights)


samplingnum = randi([1,sum(weights)]);

sample_i = 1;

if samplingnum == cumulative_weights(end)
    sample_i = numbins;
else
    while cumulative_weights(sample_i)<=samplingnum
        sample_i = sample_i+1;
    end
end


Emig = energies(sample_i);


%% calculate Efv_mean
Emig_mean = 0;
for i = 1:numbins
    Emig_mean = Emig_mean+energies(i)*weights(i);
end
Emig_mean = Emig_mean/cumulative_weights(end);


end




