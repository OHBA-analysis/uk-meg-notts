function [ statemaps ] = output_metastate_maps( hmm, metastate_hmm, FO_downsample, winsize, output_name, Ds_fnames, parcellation )

if nargin< 6
    load(hmm.filenames.concat)

    Ds_fnames=hmm.filenames.prepare;

    D=spm_eeg_load(Ds_fnames{1});
    parcellation=D.parcellation;
    
else
    %Ds_fnames=GLEAN.subspace.data;
end

logtrans=0;
embed.do=0;
embed.rectify=0;
envelope_do=0;

%%%%%
%% adjust hmm statepath to contain only metatstatepaths
hmm_meta=hmm;

if 1
    % first need to resample downsampled metatstatepaths
    ntpts=size(hmm.train.Gamma);
    ts=1/hmm.fsample:1/hmm.fsample:ntpts*1/hmm.fsample;
    ts_FO=ts(1:FO_downsample:end);

    % do winsize buffering 
    hmm_meta.statepath(:)=0;
    hmm_meta.statepath(1:winsize)=metastate_hmm.statepath(1);
    hmm_meta.statepath(end-winsize:end)=metastate_hmm.statepath(end);

    for tt=1:length(metastate_hmm.statepath)
        from=winsize+(tt-1)*FO_downsample+1;
        to=from+FO_downsample;
        hmm_meta.statepath(from:to)=metastate_hmm.statepath(tt);    
    end
    hmm_meta.K=metastate_hmm.K;
end

%%%%%


statemaps = [output_name,'metastate_map'];


clear tmp;
for subnum = 1:length(Ds_fnames)

    % compute subject's state maps
    hmm_sub = hmm_meta;
    hmm_sub.statepath = hmm_meta.statepath(hmm.subj_inds==subnum);

    % hmm_sub.train.Gamma=hmm.train.Gamma(subj_inds==subnum,:);

    %hmm_sub = rmfield(hmm_sub,'MixingMatrix');
 
    Dp = spm_eeg_load(Ds_fnames{subnum});
    
    if 1
        normalisation='voxelwise';
        output_method='pcorr';
        de_mean=1;
        diff_contrast=0;
        state_assignment='hard';
    else
        normalisation='global';
        output_method='cope';
        de_mean=0;
        diff_contrast=0;
        state_assignment='hard';
    end
    
    datap = prepare_data(Dp,normalisation,logtrans,1,embed);
    
    tmp(subnum,:,:)=osl_hmm_statemaps(hmm_sub,datap,~envelope_do,output_method,state_assignment,de_mean,diff_contrast);
            
    %statp = statp + osl_hmm_statemaps(hmm_sub,datap,~envelope_do,output_method,state_assignment,de_mean,diff_contrast);
    %statp = statp + tmp3;
end

%statp = statp ./ length(hmm.filenames.prepare);
 
statp = squeeze(mean(tmp,1));%./std(tmp,[],1));
  
statemaps=[statemaps,'_parcels'];

map = parcellation2map(statp,parcellation.file,parcellation.mask);
writenii(map,statemaps,parcellation.mask);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = prepare_data(D,normalisation,logtrans,freq_ind,embed)

% reshape trialwise data
if strcmp(D.transformtype,'TF')
    data = D(:,freq_ind,:,:);
else
    data = D(:,:,:);
end
data = reshape(data,[D.nchannels,D.nsamples*D.ntrials]);

% select only good data
good_samples = ~all(badsamples(D,':',':',':'));
good_samples = reshape(good_samples,1,D.nsamples*D.ntrials);
data = data(:,good_samples);

% log transform
if logtrans
    data = log10(data);
end

if exist('embed','var') && embed.do,    
    
    disp('Time embedding data');
    span=1/embed.centre_freq; %secs
    num_embeddings=round(span/embed.tres);        
    %lags=round(linspace(-num_embeddings/2,num_embeddings/2,num_embeddings));
    lags=round(-num_embeddings/2:num_embeddings/2);
    %lags=round(0:num_embeddings-1);

    disp(lags);
    
    [dataout,valid]=embedx(data',lags);
    dataout=dataout';
    data=randn(size(dataout,1),size(data,2))*std(squash(data(:,1:500)));
    data(:,valid)=dataout;          
       
end;

% normalisation
switch normalisation
    case 'global'
        data = demean(data,2)./std(data(:));
    case 'voxelwise'
        data = normalise(data,2);
    case 'none'
        data = demean(data,2);
end

if embed.rectify
    data=abs(data);
    
    switch normalisation
    case 'global'
        data = demean(data,2)./std(data(:));
    case 'voxelwise'
        data = normalise(data,2);
    case 'none'
        data = demean(data,2);
    end;
end

end