%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path
    
tilde='/home/mwoolrich/';
tilde='/Users/woolrich/';

cd([tilde '/homedir/scripts/osl/osl-core']);
%cd /Users/woolrich/homedir/osl/current_release/osl/osl-core

osl_startup();

addpath(genpath([tilde 'Dropbox/vols_scripts/hmm_misc_funcs']));
addpath(genpath([tilde 'Dropbox/vols_scripts/notts_ukmp/main']));

% sessionname='eo'; nohup /Applications/MATLAB_R2017b.app/bin/matlab -nodesktop -nodisplay < /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/main/main_run_hmm.m > /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/logs/main_run_hmm_$sessionname.log 2>&1 &
% sessionname='eo'; tail -f /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/logs/main_run_hmm_$sessionname.log

% more /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/main/main_run_hmm.m

% to get hmm results onto laptop:
% scp -r $M13:/Users/woolrich/homedir/vols_data/notts_ukmp/spm/sept2019_eo/hmm_1to45hz/hmm2_parc_giles_symmetric__pcdim80_voxelwise_embed13_K8_big1_dyn_modelhmm.mat /Users/woolrich/homedir/vols_data/notts_ukmp/spm/sept2019_eo/hmm_1to45hz/
% scp -r $M13:/Users/woolrich/homedir/vols_data/notts_ukmp/spm/sept2019_eo/hmm_1to45hz/hmm2_parc_giles_symmetric__pcdim80_voxelwise_embed13_K8_big1_dyn_modelhmm_output_pcorr_parcels.nii.gz /Users/woolrich/homedir/vols_data/notts_ukmp/spm/sept2019_eo/hmm_1to45hz/

% if not done previously, to get hmm prepared data onto laptop:
% scp -r $M13:/Users/woolrich/homedir/vols_data/notts_ukmp/spm/sept2019_eo/hmm_1to45hz/hmm2_parc_giles_symmetric__pcdim80_voxelwise_embed13.mat /Users/woolrich/homedir/vols_data/notts_ukmp/spm/sept2019_eo/hmm_1to45hz/
% if not done previously, to get input data into teh_groupinference_parcels 
% scp $M13:/Users/woolrich/homedir/vols_data/notts_ukmp/spm/sept2019_eo/bfnew_1to45hz/sfold_giles_symmetric_fsept2019_eo_session*.?at /Users/woolrich/homedir/vols_data/notts_ukmp/spm/sept2019_eo/bfnew_1to45hz/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

spmfilesdir='/Users/woolrich/homedir/vols_data/notts_ukmp/spm/';

session_name='eo';nsubjects=55;
%session_name='eo';nsubjects=2;
%session_name='vml';nsubjects=51;
%session_name='vml';nsubjects=2;
%session_name='vms';nsubjects=52;

preproc_name='sept2019'; 

prep_sessions_to_do=1:nsubjects;
sessions_to_exclude=[];
hmm_sessions_to_do=setdiff(prep_sessions_to_do,sessions_to_exclude);

%%

S=[];

S.do_prepare=0; 
S.do_hmm=0; 
S.do_spectral_estimation=0;
S.preproc_name=preproc_name;
S.session_name=session_name; 
S.spmfilesdir=spmfilesdir;
S.prep_sessions_to_do=prep_sessions_to_do;
S.hmm_sessions_to_do=hmm_sessions_to_do;
S.num_embeddings=13;
S.freq_range=[1 45];

S.hmm_name='';

disp(S);

[hmm, hmmfname, hmmoptions, settings_prepare, spm_files, epoched_spm_files] = run_full_hmm_pipeline(S);

disp(hmmfname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

storage_dir=[hmmfname '_store/'];

mkdir(storage_dir);

printdir=[storage_dir 'figs/'];

mkdir(printdir);

freq=[];
freq{1}=[1 45];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

osleyes([hmm.statemaps '_parcels.nii.gz']);

%%%%%%%%%%%%%%%%%%%%%%%%
%% View state maps

fsleyes([hmm.statemaps '_parcels.nii.gz']);
 
%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Quick route to plotting coherence network!

do_nnmf=1;
frontal_state=5;
states=[frontal_state];
plot_coh_graph_from_hmm(hmm,states,do_nnmf,freq,storage_dir,printdir);


%% get percentile thresholds

x=ra([hmm.statemaps '_parcels']);

clear pc
for kk=1:size(x,4)
    tmp=abs(x(:,:,:,kk));
    pc(kk)=percentile(squash(tmp,tmp),85);
end
pc

%%

fout = osl_resample_nii([hmm.statemaps '_parcels'], [hmm.statemaps '_parcels_2mm'], 2, 'nearestneighbour', [osldir '/std_masks/MNI152_T1_2mm_brain']);

fslview(fout);

workbenchdir=['/Applications/workbench/bin_macosx64'];
osl_render4D([fout '.nii.gz'],printdir,workbenchdir,'nearestneighbour',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot temporal statistics

figure
stats = osl_hmm_stats(hmm,'do_plots');

print([printdir 'temporal_stats'],'-dpng');

%

ts=0:0.1:12
clear lts
figure;
for ii=1:length(stats),
    lts(ii,:)=hist(stats(ii).LifeTimes,ts);
    ho;
    plot(ts,lts(ii,:)/sum(lts(ii,:)),'LineWidth',2);
end;
plot4paper('time (s)','');
title('Life times');
print([printdir 'lifetimes'],'-dpng');

ts=0:0.02:12
clear lts 
figure;
for ii=1:length(stats),
    lts(ii,:)=hist(stats(ii).Intervals,ts);
    ho;
    plot(ts,lts(ii,:)/sum(lts(ii,:)),'LineWidth',2);
   
end;
plot4paper('time (s)','');
title('Interval Times');
print([printdir 'intervaltimes'],'-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

clear lts mils longones legs
figure;
for ii=1:length(stats),
    lts{ii}=stats(ii).LifeTimes;
    mils{ii}=stats(ii).Intervals;  
    
    for tt=1:10
        longones(ii,tt)=mean(mils{ii}>tt)*100;
    end
    legs{ii}=['S' num2str(ii)];
    plot(longones(ii,:),get_cols(ii),'LineWidth',3);
    ho;
end
legend(legs);
plot4paper('Interval time, t (s)','Percentage of intervals longer than time, t');
print([printdir 'intervaltimes_gt_time'],'-dpng');

%

figure;
[h]=distributionPlot(lts,'histOpt',2,'color','c');
h{2}(1).Marker='.';
h{2}(1).MarkerSize=40;
h{2}(1).Color='r';
h{2}(2).Marker='.';
h{2}(2).MarkerSize=40;
h{2}(2).Color='k';

%[h,L,~,~,bw]=violin(lts); L.Visible='off';
plot4paper('State #','Dwell Time (s)');
a=axis;
axis([a(1) a(2) 0 0.3]);
fig=gcf;
fig
print([printdir 'violin_lifetimes'],'-dpng');

figure;
[h]=distributionPlot(mils,'histOpt',2,'color','c');
h{2}(1).Marker='.';
h{2}(1).MarkerSize=40;
h{2}(1).Color='r';
h{2}(2).Marker='.';
h{2}(2).MarkerSize=40; 
h{2}(2).Color='k';
%[h,L]=violin(mils); L.Visible='off';
plot4paper('State #','Interval Time (s)');
a=axis;
axis([a(1) a(2) 0 2]);
fig=gcf;
fig
print([printdir 'violin_intervaltime'],'-dpng');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot state evoked responses

[gamma_epoched, Ds]=teh_epoch_statetcs(hmm, epoched_spm_files);

%

clear gamma_epoched_conds
for subnum=1:length(hmm.data_files)
    for cc =1:length(Ds{1}.condlist)
        gamma_epoched_conds(cc,subnum,:,:)=mean(Ds{subnum}(find(strcmp(Ds{subnum}.chantype,'CLASSPR')),:,Ds{subnum}.indtrial(Ds{subnum}.condlist{cc},'good')),3);
    end
end

%%

%%for cc =1:length(Ds{1}.condlist)

xrange=[-5 5];
cc=1;baseline_time=[-5 -2];
cc=2;baseline_time=[-5 -2];

fo=squeeze(mean(gamma_epoched_conds(cc,:,:,:),2));

winsize=1;
legs=[];
for kk=1:size(fo,1),
    fo(kk,:)=moving(fo(kk,:),winsize);
    legs{kk}=num2str(kk);
end
figure;plot(Ds{1}.time,fo,'LineWidth',2);
title(Ds{1}.condlist{cc})
a=axis; a(1:2)=xrange; axis(a);
legend(legs);
print([printdir 'fo_epoched_cond' num2str(cc)],'-dpng');

baseline_inds=find(Ds{1}.time>baseline_time(1) & Ds{1}.time<baseline_time(2));

legs=[];
clear fo_bc fo_bcs
for kk=1:size(fo,1),
    fo_bc=fo(kk,:);
    fo_bc=fo_bc-mean(fo_bc(baseline_inds));
    fo_bc=moving(fo_bc,winsize);
    fo_bcs(kk,:)=fo_bc;
    legs{kk}=num2str(kk);
end
figure;plot(Ds{1}.time,fo_bcs,'LineWidth',2);
a=axis; a(1:2)=xrange; axis(a);
title(Ds{1}.condlist{cc})
legend(legs);
print([printdir 'fo_epoched_bc_cond' num2str(cc)],'-dpng');

figure;
for kk=1:hmm.K
    subplot(2,ceil(hmm.K/2),kk)
    plot(Ds{1}.time,fo(kk,:),'LineWidth',2);
    title([Ds{1}.condlist{cc} ', k=' num2str(kk)]);
    a=axis; a(1:2)=xrange; axis(a);
end

print([printdir 'fo_epoched_sep_cond' num2str(cc)],'-dpng');

%%

for kk=1:hmm.K
    figure;
    
    plot(Ds{1}.time,fo(kk,:),'LineWidth',2);
    title([Ds{1}.condlist{cc} ', k=' num2str(kk)]);

    a=axis; a(1:2)=xrange; axis(a);
    print([printdir 'fo_epoched_cond' num2str(cc) '_k' num2str(kk)],'-dpng');
end



%%%%%%%%%%%%%%%%%%%%%%%%
%% superstates?

clear cnt
for ss=1:length(hmm.epoched_statepath_sub)
   tmp= hmm.epoched_statepath_sub{ss};
   for kk=1:hmm.K
    cnt(ss,kk)=sum(tmp==kk);
   end
   cnt(ss,:)=cnt(ss,:)/length(tmp);
end
figure; imagesc(corrcoef(cnt));
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%
%% trans probs mat

tmp=hmm.P;
tmp=tmp-diag(diag(tmp));
figure;imagesc(tmp);
print([printdir 'transprob'],'-dpng');

tmp=hmm.P;
%tmp=tmp-diag(diag(tmp));
Z=linkage(tmp','ward')
figure;dendrogram(Z)
print([printdir 'dendogram'],'-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot continuous statepath
figure
option=[];
option.mode='separate';
option.win=[];

hmmtmp=rmfield(hmm,'epoched_statepath_sub');
options=[]; options.time_range=[1 1000];
osl_hmm_plotstatepath(hmmtmp,options);
set(gcf,'Position', [63, 424, 1367, 374])

fig = gcf;
fig.PaperPositionMode = 'auto';
print([printdir 'statepath_cont'],'-dpng');

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
    load([storage_dir '/state_netmats_mt_' num2str(floor(S.netmat_method_options.reg)) ...
    '_vn' num2str(S.netmat_method_options.var_normalise) '_' S.assignment  '_' ...
    'global' num2str(S.global_only)],'state_netmats_mt');

end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate spectra and cross spectra using multitaper on each session separately

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot node and edge spectral properties

node_pairs={};
node_pairs{1}=[17 18];
node_pairs{1}=[37 38];


if 0
    
    node_pairs{1}=[9 10];
    node_pairs{2}=[19 20];
    node_pairs{3}=[9 19];
    node_pairs{4}=[37 38];
    node_pairs{5}=[17 18];
    node_pairs{6}=[38 17];
    node_pairs{7}=[38 18];

end

print_plot=1;
for nn=1:length(node_pairs)
       
    dodiff=0;
    meas='coh';
        
    legs=hmm_state_plot_coh_concat(state_netmats_mt, node_pairs{nn}, meas, dodiff, print_plot);
    
    state_to_plot_all_subjects_for=1;
    hmm_state_plot_coh(state_netmats_mtsess, node_pairs{nn}, meas, state_to_plot_all_subjects_for);
    
    subplot(2,2,1);title(['reg' num2str(floor(S.netmat_method_options.reg))  ', vn' num2str(S.netmat_method_options.var_normalise)]);
    print([printdir 'pairwise_spectra_' num2str(node_pairs{nn}(1)) '_' num2str(node_pairs{nn}(2))],'-dpng');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dodgey state?

clear numpts;
for ss=1:length(state_netmats_mtsess)
    for kk=1:hmm.K
        numpts(ss,kk)=state_netmats_mtsess{ss}.state{kk}.ntpts;
    end
end

dodgey_state=6;
figure;plot(numpts(:,dodgey_state));

%%

dodgey_state=6;

for ss=1:length(state_netmats_mtsess)
    state_netmats_mtsess{ss}.state{dodgey_state}=state_netmats_mtsess{ss}.state{1};
end


state_netmats_mt{1}.state{dodgey_state}=state_netmats_mt{1}.state{1};

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
        fres=2;
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

do_run=1;
S=[];
S.psds=psds_fbs;

% for alpha, beta, delta
S.maxP=4; S.maxPcoh=4; 

% for wideband
S.maxP=1; S.maxPcoh=2; 
 
if do_run
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
%% plot nnmf psd spectra
clear a statset
for pp=1:S.maxP
    figure; 
    plot(freqBandmn, squeeze(spectral_nnmf_res.nnmf_psd_specs(pp,:)), get_cols(pp), 'LineWidth',4);
    plot4paper('freq (Hz)', '');
    set(gca,'YtickLabels',[]);

    print([printdir 'nnmf' num2str(pp) '_specs'],'-dpng');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write nnmf psd maps by freq for view in fsleyes

clear node_maps_fname fout;
for pp=1:S.maxP
    
    nodemaps=squeeze(spectral_nnmf_res.nnmf_psd_maps(:,pp,:))';
    node_maps_fname{pp}=[printdir 'fac' num2str(pp) '_map'];
    sres=nii.get_spatial_res(hmm.parcellation.mask);
    map = parcellation2map(nodemaps,hmm.parcellation.file,sres(1));
    nii.quicksave(map,node_maps_fname{pp},sres(1));
    
    fout_lowres{pp}=[node_maps_fname{pp} '.nii.gz'];
    
    fout{pp} = nii.resample(fout_lowres{pp}, [node_maps_fname{pp} '_2mm'], 2,'interptype', 'nearest');
    fout{pp}=[fout{pp} '.nii.gz'];
end

% num of fout files is num of factors in nnmf
% each fout is 4D nii with 4th dim being states
fsleyes(fout); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% open nnmf psd maps by freq in workbench

for pp=1:maxP-1
    workbenchdir=['/Applications/workbench/bin_macosx64'];
    osl_render4D([fout{pp} '.nii.gz'],printdir,workbenchdir,'nearestneighbour',1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write nnmf psd maps by state for nnmf

clear node_maps_fname fout;
cmd=['fslview ' osldir '/std_masks/MNI152_T1_8mm_brain.nii.gz '];
for kk=1:NK

    nodemaps=squeeze(spectral_nnmf_res.nnmf_psd_maps(kk,:,:))';
    node_maps_fname{kk}=[printdir 'k' num2str(kk) '_map'];
    sres=nii.get_spatial_res(hmm.parcellation.mask);
    map = parcellation2map(nodemaps,hmm.parcellation.file,sres(1));
    
    nii.quicksave(map,node_maps_fname{kk},sres(1));
    
    fout_lowres{kk}=[node_maps_fname{kk} '.nii.gz'];
    
    fout{kk} = nii.resample(fout_lowres{kk}, [node_maps_fname{kk} '_2mm'], 2,'interptype', 'nearest');
    fout{kk}=[fout{kk} '.nii.gz'];
end

% num of fout files is num of factors in nnmf
% each fout is 4D nii with 4th dim being states
fsleyes(fout); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% open nnmf psd maps by state in workbench

for kk=1:NK
    nodemaps=squeeze(spectral_nnmf_res.nnmf_psd_maps(kk,:,:))';
    kk
    percentile(squash(nodemaps),85)
    
    workbenchdir=['/Applications/workbench/bin_macosx64'];
    osl_render4D([fout{kk} '.nii.gz'],printdir,workbenchdir,'nearestneighbour',1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot nnmf coh spectra
figure;
for pp=1:S.maxPcoh
    plot(freqBandmn, squeeze(spectral_nnmf_res.nnmf_coh_specs(pp,:)), get_cols(pp), 'LineWidth',2);
    ho;
end
plot4paper('freq (Hz)', 'Power (a.u.)');

print([printdir 'nnmf_coh_maxp' num2str(S.maxPcoh) '_specs'],'-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at distribution over all states and edges using mixture model thresholding
clear data
for kk=[1:NK]
    
    for pp=1:S.maxPcoh-1
        
        graph=abs(squeeze(spectral_nnmf_res.nnmf_coh_maps(kk,pp,:,:)));
        tmp=squash(triu(graph));
        inds2=find(tmp~=0);
        data(kk,pp,:)=tmp(inds2);       
    end
    
end

S2=[];
S2.data=squash(data);
S2.do_fischer_xform=true;
graph_ggm=teh_graph_gmm_fit(S2);
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% produce nnmf coh brain vids using mixture model thresholding

for kk=[7]
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
        S2.do_fischer_xform=true;
        S2.pvalue_th=0.025;
        graph_ggm=teh_graph_gmm_fit(S2);
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print([printdir 'graph_gmm_fit_k' num2str(kk) '_p' num2str(pp) '_maxp' num2str(S.maxPcoh) '_p' num2str(pp) ],'-dpng');


        if 1
            th=graph_ggm.normalised_th;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maps as nii

for pp=1:S.maxPcoh-1
    
    clear nodemaps;

    for kk=1:NK, 
        graph=squeeze(sum(abs(spectral_nnmf_res.nnmf_coh_maps(kk,pp,:,:)),4));
        graph=normalise(graph);
        nodemaps(:,kk)=graph;     
    end

    node_maps_fname=[printdir 'nodeweight_map_nnmf' num2str(pp)];
    map = parcellation2map(nodemaps,hmm.parcellation.file,sres(1));
    nii.quicksave(map,node_maps_fname,sres(1));    

    %fslview(node_maps_fname);

    fout = nii.resample([node_maps_fname '.nii.gz'], [node_maps_fname '_2mm.nii.gz'], 2,  'interptype', 'nearest');
    
    %fslview(fout);

    workbenchdir=['/Applications/workbench/bin_macosx64'];
    osl_render4D([fout '.nii.gz'],printdir,workbenchdir,'nearestneighbour',1);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% amp vs connectivity using subject-wise coherence

%% need to apply precomputed nnmf specs to subject specific psds to get subject specific coh maps and amp maps

do_group_av=1;

% plot precomputed spec
figure;
for pp=1:S.maxPcoh   
    plot(spectral_nnmf_res.nnmf_coh_specs(pp,:),get_cols(pp),'Linewidth',2);ho;
end
legend(get_cols)

if do_group_av
    % do group average
    num_su=1;
    maps_coh_subj=zeros(num_su,NK,size(spectral_nnmf_res.nnmf_coh_specs,1),num_nodes,num_nodes);
    maps_auto_subj=zeros(num_su,NK,size(spectral_nnmf_res.nnmf_coh_specs,1),num_nodes);
    
    S2=[];
    S2.psds=psds_fbs;
    S2.maxP=size(spectral_nnmf_res.nnmf_coh_specs,1); 
    S2.maxPcoh=S2.maxP;
    S2.fixed_psd_specs=spectral_nnmf_res.nnmf_coh_specs;
    S2.fixed_coh_specs=spectral_nnmf_res.nnmf_coh_specs;
    S2.do_plots=false;
    spectral_nnmf_res_sub=teh_spectral_nnmf(S2);

    maps_coh_subj(1,:,:,:,:)=spectral_nnmf_res_sub.nnmf_coh_maps(1:NK,:,:,:);    
    maps_auto_subj(1,:,:,:)=spectral_nnmf_res_sub.nnmf_psd_maps(1:NK,:,:);

    
else    
    num_su=size(psds_fbs,1);
    maps_coh_subj=zeros(num_su,NK,size(spectral_nnmf_res.nnmf_coh_specs,1),num_nodes,num_nodes);
    maps_auto_subj=zeros(num_su,NK,size(spectral_nnmf_res.nnmf_coh_specs,1),num_nodes);

    for ss=1:num_su,
        S2=[];
        S2.psds=psds_fbs(ss,:,:,:,:);
        S2.maxP=size(spectral_nnmf_res.nnmf_coh_specs,1); 
        S2.maxPcoh=S2.maxP;
        S2.fixed_psd_specs=spectral_nnmf_res.nnmf_coh_specs;
        S2.fixed_coh_specs=spectral_nnmf_res.nnmf_coh_specs;
        S2.do_plots=false;
        spectral_nnmf_res_sub=teh_spectral_nnmf(S2);

        maps_coh_subj(ss,:,:,:,:)=spectral_nnmf_res_sub.nnmf_coh_maps(1:NK,:,:,:);    
        maps_auto_subj(ss,:,:,:)=spectral_nnmf_res_sub.nnmf_psd_maps(1:NK,:,:);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots

if 1
    % show over subjs
    nodeweights=(mean(mean(abs(maps_coh_subj),5),4));    
    amp_persubj=(mean((maps_auto_subj),4));
else
    % show over nodes
    nodeweights=permute((mean(mean(abs(maps_coh_subj),5),1)),[3 1 2]);
    amp_persubj=permute((mean((maps_auto_subj),1)),[3 1 2]);
end

neworder=1:NK;
[~,neworder]=sort(mean(amp_persubj,1));
    
for pp=1:S.maxPcoh-1

    figure;
    cols={'ro','go','bo','mo','co','yo','rd','gd','bd','md','cd','yd'};

    for kk=1:NK
        plot(amp_persubj(:,neworder(kk),pp),nodeweights(:,neworder(kk),pp),cols{kk},'MarkerSize',6);
        legs{kk}=['S' num2str(neworder(kk))];
        ho;
    end
    plot4paper('Amplitude (AU)','Node weights (AU)');
    legend(legs);
    print([printdir 'power_vs_coh_nodeweight_bystate_maxp' num2str(S.maxPcoh) '_p' num2str(pp) '_av' num2str(do_group_av)],'-dpng');
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get cleaned coherences
% needs spectral_nnmf_res to contain to 2 modes, the 2nd is assumed to be
% noise

coh_clean=squeeze(mean(spectral_nnmf_res.coh,1));

nROIs=size(coh_clean,3);
nF=size(coh_clean,2);

for kk=1:NK
    cohtmp=squeeze(mean(spectral_nnmf_res.coh(:,kk,:,:,:),1)); 
    cohtmp=reshape(cohtmp,[nF,nROIs*nROIs]);
    
    % regress noise mode spectra out of cohtmp
    noise_ev=spectral_nnmf_res.nnmf_coh_specs(2,:)';

    %for jj=1:size(cohtmp,2)
        %beta=lsqnonneg(?,cohtmp(:,jj));
        %cohtmp(:,jj)=cohtmp(:,jj)-noise_ev*beta;
    %end
    
    coh_clean(kk,:,:,:)=reshape(cohtmp,[nF,nROIs,nROIs]);
end

figure;
plot(squeeze(coh_clean(:,:,9,10))');
figure;
plot(squeeze(mean(spectral_nnmf_res.coh(:,:,:,37,38),1))');
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now do it with prespecifed freq bands rather than using NNMF

freq_range=[];
freq_range{1}=[1 30];

%freq_range{1}=[1 7];
%freq_range{2}=[7 13];
%freq_range{3}=[13 20];

%%

cmd=['fslview ' OSLDIR '/std_masks/MNI152_T1_8mm_brain.nii.gz '];

clear node_maps_fname fout;
for ff=1:length(freq_range)

    freq_inds=find(freqBandmn>freq_range{ff}(1) & freqBandmn<freq_range{ff}(2));

    %%%
    % DO PSD
    clear psd;
    for kk=1:NK
        tmp=squeeze(mean(abs(auto_spectra_comps(1,kk,freq_inds,:)),3));
        psd(kk,:)=(tmp-mean(tmp))./std(tmp);
    end
    
    % write maps
    nodemaps=psd';
    node_maps_fname{ff}=[printdir 'freq' num2str(ff) '_map'];
    map = parcellation2map(nodemaps,hmm.parcellation.file,hmm.parcellation.mask);
    writenii(map,node_maps_fname{ff},hmm.parcellation.mask);
    fout{ff} = osl_resample_nii(node_maps_fname{ff}, [node_maps_fname{ff} '_2mm'], 2, 'nearestneighbour', [osldir '/std_masks/MNI152_T1_2mm_brain']);
    cmd=[cmd node_maps_fname{ff} ' '];

end
cmd=[cmd ' &'];

    %fslview(fout);

%% maps for psd from freq bands
% open workbench
ff=1;
workbenchdir=['/Applications/workbench/bin_macosx64'];
osl_render4D([fout{ff} '.nii.gz'],printdir,workbenchdir,'nearestneighbour',1)

%runcmd(cmd);

%% graphs for coh from freq bands

for ff=1:length(freq_range)

    freq_inds=find(freqBandmn>freq_range{ff}(1) & freqBandmn<freq_range{ff}(2));
    
    %%%
    % DO COH
    clear graphs;
    for kk=1:NK
        psd=squeeze(mean(abs(coh_comps(1,kk,freq_inds,:,:)),3));
        psd=psd-diag(diag(psd));
        graphs(kk,:,:)=abs(psd);
    end

    % threshold network matrix 
    for kk=[4]

        graph=squeeze(graphs(kk,:,:));
        tmp=squash(triu(graph));
        inds=find(tmp~=0);
        tmp=tmp(inds);
        tmp2=log((tmp)./(1-tmp));
        tmp3=(tmp2-mean(squash(tmp2)))./std(squash(tmp2));
        th=2;
        graph=tmp3;
        graph(graph<th)=NaN;
        graphmat=nan(num_nodes, num_nodes);
        graphmat(inds)=graph;
        graph=graphmat;
        
        % get node co-ordinates parcellation
        parc_file=[tilde '/Dropbox/vols_scripts/parcellations/fmri_d100_parcellation_with_PCC_reduced_2mm'];
        options.parcellation.file = [parc_file '_ss5mm_ds8mm'];   
        parcelFile = options.parcellation.file;
        spatialRes = 8; 
        spatialMap = nii.quickread(parcelFile, spatialRes); 
        mni_coords = find_ROI_centres(spatialMap, spatialRes, 0); 

        % plot 
        nROIs = size(graph,2);
        colorLims = [2 4]; 
        sphereCols = repmat([30 144 255]/255, nROIs, 1); 
        edgeLims = [4 8];
        figure;%('Color', 'w', 'Name', [figureStem ' ' BAND_NAMES{iBand} ' band']);
        osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], 0, sphereCols, edgeLims); 
        %

        fileName=[printdir 'brain_vid_' num2str(freq_range{ff}(1)) 'hz_to_' num2str(freq_range{ff}(2)) 'hz_k' num2str(kk)];

        make_brain_vid(fileName)

    end
end






