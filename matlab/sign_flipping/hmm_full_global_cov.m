function [netmats] = hmm_full_global_cov(Ds, S)

logtrans = 0;

normalisation    = S.concat.normalisation;
embed            = S.concat.embed;

nsubs = length(Ds);

for subnum = 1:nsubs
    disp(['Preparing subject ' num2str(subnum)]);
    Dp = spm_eeg_load(Ds{subnum});
    embed.tres=1/Dp.fsample;
            
    % returns data as num_nodes x num_embeddings x ntpts
    [datap, num_embeddings] = prepare_data(Dp, normalisation, logtrans, embed);
       
    netmats{subnum}.global = netmat_cov(datap, S);
    netmats{subnum}.num_embeddings = num_embeddings;
end    
