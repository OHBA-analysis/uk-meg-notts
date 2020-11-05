function plot_coh_graph_from_hmm(hmm,states,do_nnmf,freq,storage_dir,printdir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate spectra and cross spectra using multitaper on concatenated data

do_run=0;

S=[];
S.parcellated_filenames=hmm.data_files;
S.normalisation='voxelwise';
S.assignment='hard';
S.global_only=false;
S.embed.do=0;
S.embed.rectify=false;

S.netmat_method=@netmat_spectramt;
S.netmat_method_options.fsample=hmm.fsample;
S.netmat_method_options.fband=freq{1};
S.netmat_method_options.type='coh';
S.netmat_method_options.full_type='full';
S.netmat_method_options.var_normalise=false;
S.netmat_method_options.reg=2; % higher is less reg
S.netmat_method_options.order=0;

if do_run

    [ state_netmats_mt ] = hmm_state_netmats_teh_concat( hmm, S );

    save([storage_dir '/state_netmats_mt' num2str(floor(S.netmat_method_options.reg)) ...
        '_vn' num2str(S.netmat_method_options.var_normalise) '_' S.assignment '_' ...
        'global' num2str(S.global_only)], '-v7.3', 'state_netmats_mt');
else
    load([storage_dir '/state_netmats_mt' num2str(floor(S.netmat_method_options.reg)) ...
    '_vn' num2str(S.netmat_method_options.var_normalise) '_' S.assignment  '_' ...
    'global' num2str(S.global_only)],'state_netmats_mt');

end

%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate spectra and cross spectra using multitaper on each session separately

do_run=0;

S=[];
S.parcellated_filenames=hmm.data_files;
S.normalisation='voxelwise';
S.assignment='hard';
S.global_only=false;
S.embed.do=0;
S.embed.rectify=false;

S.netmat_method=@netmat_spectramt;
S.netmat_method_options.fsample=hmm.fsample;
S.netmat_method_options.fband=freq{1};
S.netmat_method_options.type='coh';
S.netmat_method_options.full_type='full';
S.netmat_method_options.var_normalise=false;
S.netmat_method_options.reg=2;
S.netmat_method_options.order=0;

if do_run

    [ state_netmats_mtsess ] = hmm_state_netmats_teh( hmm, S );

    save([storage_dir '/state_netmats_mtsess_' num2str(floor(S.netmat_method_options.reg)) ...
        '_vn' num2str(S.netmat_method_options.var_normalise) '_' S.assignment '_' ...
        'global' num2str(S.global_only)], '-v7.3', 'state_netmats_mtsess');
    
else
    
    load([storage_dir '/state_netmats_mtsess_' num2str(floor(S.netmat_method_options.reg)) ...
    '_vn' num2str(S.netmat_method_options.var_normalise) '_' S.assignment  '_' ...
    'global' num2str(S.global_only)],'state_netmats_mtsess');

end

disp('coh loaded');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spatial map FC by freq 

use_sess=1; % NNMF works much better on concatenated data over subjects (rather than on group averaged spectra)

if ~use_sess
    state_netmats=state_netmats_mt;
else
    state_netmats=state_netmats_mtsess;
end
    
NK=length(state_netmats{1}.state);
num_nodes=size(state_netmats{1}.state{1}.netmat,1);
num_freqs=length(state_netmats{1}.state{1}.spectramt.f);
nsubjects=length(state_netmats);

psds=zeros(nsubjects,NK+1,num_freqs,num_nodes,num_nodes);
for ss=1:length(state_netmats)
    for kk=1:length(state_netmats{1}.state)           
        psds(ss,kk,:,:,:)=state_netmats{ss}.state{kk}.spectramt.psd;                               
    end
    % add global on the end
    psds(ss,kk+1,:,:,:)=state_netmats{ss}.global.spectramt.psd;
end


% compute freq bands
fs=state_netmats{1}.state{1}.spectramt.f;
do_freq_downsample=true;

freqBands=[];
freqBandmn=[];
if 0
    freqBands{1}=[1 4];
    freqBands{2}=[4 8];
    freqBands{3}=[8 12];
    freqBands{4}=[12 20];
    freqBands{5}=[20 30];
    freqBands{6}=[30 45];
else
    if do_freq_downsample
        fres=1;
        from = 1;
        to = from + fres;  
        ff=1;
        while to<46                
            freqBands{ff}=[from to]; 
            freqBandmn(ff)=mean(freqBands{ff});
            ff=ff+1;
            from = to;
            to = from + fres;        
        end
    else
        clear freqBands;
        freqBandmn=state_netmats{1}.state{1}.spectramt.f;
    end
end

nfreqbins=length(freqBandmn);

if do_freq_downsample
    clear psds_fbs;
    for ss=1:length(state_netmats)
        for kk=1:NK+1, 
            for ff=1:length(freqBands),

                finds=find(fs>freqBands{ff}(1) & fs < freqBands{ff}(2));

                psds_fbs(ss,kk,ff,:,:) = squeeze(mean(psds(ss,kk,finds,:,:),3)); 

            end    
        end
    end
else
    psds_fbs=psds;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run or load NNMF

clear spectral_nnmf_res;

S=[];
S.psds=psds_fbs;

% for alpha, beta, delta
S.maxP=4; S.maxPcoh=4; 

% for wideband
S.maxP=1; S.maxPcoh=2; 
 
if do_nnmf
    if exist('spectral_nnmf_res','var')
        if isfield(spectral_nnmf_res,'nnmf_psd_specs')
            S.nnmf_psd_specs=spectral_nnmf_res.nnmf_psd_specs;
            S.nnmf_psd_maps=spectral_nnmf_res.nnmf_psd_maps;
        end    
    end

    spectral_nnmf_res=teh_spectral_nnmf(S);
    
    save([storage_dir '/spectral_nnmf_res_maxP' num2str(S.maxP) '_bysess' num2str(use_sess)], '-v7.3', 'spectral_nnmf_res');
else
    load([storage_dir '/spectral_nnmf_res_maxP' num2str(S.maxP) '_bysess' num2str(use_sess)], 'spectral_nnmf_res');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% produce nnmf coh brain vids using mixture model thresholding

for kk=states
%for kk=[1,4,7,8]%[1:NK]
%for kk=[1:NK]
     
    kk
    
    for pp=1:S.maxPcoh-1
        
        graph=abs(squeeze(spectral_nnmf_res.nnmf_coh_maps(kk,pp,:,:)));
        tmp=squash(triu(graph));
        inds2=find(tmp>1e-10);
        data=tmp(inds2);
        
        S2=[];
        S2.data=squash(data);
        S2.do_fischer_xform=false;
        S2.pvalue_th=0.005; 
        graph_ggm=teh_graph_gmm_fit(S2);
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print([printdir 'graph_gmm_fit_k' num2str(kk) '_p' num2str(pp) '_maxp' num2str(S.maxPcoh) '_p' num2str(pp) ],'-dpng');

        if 1
            th=graph_ggm.normalised_th
            graph=graph_ggm.data';

            graph(graph<th)=NaN;
            graphmat=nan(num_nodes, num_nodes);
            graphmat(inds2)=graph;
            graph=graphmat;

            % get node co-ordinates parcellation
            parcelFile = hmm.parcellation.file;
            spatialRes = 8; 
            spatialMap = nii.quickread(parcelFile, spatialRes); 
            mni_coords = find_ROI_centres(spatialMap, spatialRes, 0); 

            % plot 
            nROIs = size(graph,2);
            colorLims = [th th+1]; 
            sphereCols = repmat([30 144 255]/255, nROIs, 1); 
            edgeLims = [4 8];

            figure;%('Color', 'w', 'Name', [figureStem ' ' BAND_NAMES{iBand} ' band']);
            osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], 0, sphereCols, edgeLims); 


            fileName=[printdir 'brain_vid_maxp' num2str(S.maxP) '_nnmf' num2str(pp) '_k' num2str(kk)];

           % make_brain_vid(fileName)
        end
    end
end   
