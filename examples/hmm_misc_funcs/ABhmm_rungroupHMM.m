function [HMMresults,statemaps,statePath,sub_statePath] = ABhmm_rungroupHMM(BFfiles,hmmdir,todo,env_settings,hmm_settings)
% Runs a group HMM analysis on amplitude envelopes of source reconstructed
% MEG data using the OAT beamformer
% ABhmm_rungroupHMM(BFfiles,hmmdir,todo,env_settings,hmm_settings)
%
% INPUTS:
%
% BFfiles - list of OAT source recon results (oat_stage_results)
% 
% hmmdir  - directory within which to save HMM results
%
% todo    - binary flags for each stage of the analysis [0/1]:
%           .envelope - computes amplitude envelopes from source
%                        reconstructed data
%           .prepare  - concatenates data and performs PCA dimensionality
%                        reduction
%           .infer    - infers the HMM
%
% env_settings - settings for envelope computation
%                .windowsize - window size for moving average filter [s]
%                               (default 0.1)
%                .overlap    - overlap for moving average filter [0-1]
%                               (default 0.75)
%                .norm_subs  - apply subject variance normalisation [0/1]
%                               (default 1)
%             
%
% hmm_settings - settings for HMM inference
%                .pcadim  - dimensionality to use (default 40)
%                .nstates - number of states to infer (default 8)
%                .nreps   - number of repeat inferences to run (default 5)
%
% AB 2013

global OSLDIR

HMMresults = [];
statemaps  = [];

BFnames = cell(size(BFfiles));
for f = 1:numel(BFfiles)
  [~,BFnames{f},~] = fileparts(BFfiles{f});
end

if ~isdir(hmmdir); mkdir(hmmdir); end


% Default todo settings
try todo.envelope = todo.envelope; catch, todo.envelope = 1; end
try todo.prepare  = todo.prepare;  catch, todo.prepare  = 1; end
try todo.infer    = todo.infer;    catch, todo.infer    = 1; end

% Default envelope settings
try windowsize = env_settings.windowsize; catch, windowsize = 0.1;  end
try overlap    = env_settings.overlap;    catch, overlap    = 0.75; end
try norm_subs  = env_settings.norm_subs;  catch, norm_subs  = 1;    end

% Default HMM settings
try  pcadim  = hmm_settings.pcadim;  catch, pcadim  = 40; end
try  nstates = hmm_settings.nstates; catch, nstates = 8; end
try  nreps   = hmm_settings.nreps;   catch, nreps   = 5;  end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    C O M P U T E   E N V E L O P E S                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.envelope
  
  for subnum = 1:length(BFfiles)
    
    [pathstr,filestr] = fileparts(BFfiles{subnum});
  
    disp(['Computing envelope data for ' filestr]);
    
    load(BFfiles{subnum});   
    oat_stage_results.source_recon.dirname = pathstr;
    
    tbad = ~oat_stage_results.samples2use;
    
    recon = AB_get_source_timecourses(oat_stage_results,'noreject','norecon');
    
    % Apply hilbert transform in sensor space:
    
    %%%%%%%%%%%%%
    %% MWW addition
    ntrials=size(recon.sensor_data,3);
    if ntrials>1
        %recon.sensor_data=reshape(recon.sensor_data,[size(recon.sensor_data,1), size(recon.sensor_data,2)*size(recon.sensor_data,3)]);
        %recon.sensor_data = transpose(hilbert(recon.sensor_data'));
        for tri=1:ntrials,
            recon.sensor_data(:,:,tri)=transpose(hilbert(recon.sensor_data(:,:,tri)'));
        end;
    else
        recon.sensor_data = transpose(hilbert(recon.sensor_data'));    
    end;
    %%%%%%%%%%%%%%
    
    % Preallocation
    Nvoxels  = length(recon.weights);
    winsize  = windowsize * oat_stage_results.BF.data.D.fsample;
    t        = oat_stage_results.BF.data.D.time;
    t(tbad)  = nan;
    
    % Apply temporal smoothing and downsampling - multiple orientations will be concatenated
    if winsize ~= 0
      env = zeros(Nvoxels, length(osl_movavg(t,t,winsize,overlap,0)), ntrials);
    else
      env = zeros(Nvoxels, length(t), ntrials);
    end
    
    ft_progress('init', 'etf');
    for vox=1:Nvoxels
      ft_progress(vox/Nvoxels);
      
      tdata = AB_get_source_timecourses_recon(recon,vox);
      tdata = abs(tdata{1});           
      
      if winsize ~= 0
        for tri=1:ntrials,
            HE = sqrt(sum(tdata.^2,1)); % Vector beamformer compatible
            HE(:,tbad) = nan;
 
            [env(vox,:,tri),t_avg] = osl_movavg(permute(HE(1,:,tri),[2 1 3])',t,winsize,overlap,0);
        end;
      else
        HE = sqrt(sum(tdata.^2,1)); % Vector beamformer compatible
        HE(:,tbad) = nan;
 
        env(vox,:,:) = HE;
        t_avg = t;
      end
      
    end
    ft_progress('close')
    
    t=t_avg;
    
    clear recon
    env(:,isnan(env(1,:))) = [];
    
    for tri=1:size(env,3),
        env2(:,isnan(env(1,:,tri)),tri) = [];
    end;
    
    savefile = fullfile(hmmdir,[BFnames{subnum} '_movavg' num2str(windowsize*1000) 'ms']);
    save(savefile,'env','t','-v7.3')
      
    clear oat_stage_results
    dos(['rm ' hmmdir '*.nii.gz'])
    
  end
  
end % todo.envelope




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    P R E P A R E   H M M   D A T A                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(BFfiles) == 1
  savefile_hmm = [hmmdir,BFnames{1} '_movavg' num2str(windowsize*1000) 'ms' '_hmm'];
else
  savefile_hmm = [hmmdir,'multisub' '_movavg' num2str(windowsize*1000) 'ms' '_hmm'];
end


if todo.prepare || (todo.infer && ~exist([savefile_hmm 'data.mat'],'file'))
  
  % Load subjects and concatenate:
  env_concat = [];
  subj_inds = [];
  C = 0;
  for subnum = 1:length(BFfiles)
    load([hmmdir BFnames{subnum} '_movavg' num2str(windowsize*1000) 'ms']);
    
    if size(env,3)>1, % More than 1 trial so concatenate them, MWW
        env_reshaped=reshape(env,[size(env,1) size(env,2)*size(env,3)]);
    end;
    
    if norm_subs
      env_reshaped = demean(env_reshaped,2)./env_reshaped(env(:));
    end
    
    env_concat = [env_concat,env_reshaped];
    subj_inds = [subj_inds,subnum*ones(1,size(env_reshaped,2))];
    
    C = C + env_reshaped*permute(env_reshaped,[2,1]);
    
  end
  C = C ./ (length(subj_inds)-1);
  clear env_reshaped env
  
  
  % PCA + whitening - below is equivalent to fastica whitening code but much faster
  [allsvd,Apca] = eigdec(C,pcadim); % pca([tpts x voxels],pcadim)
  whiteningMatrix = diag(1./sqrt(allsvd)) * Apca';
  hmmdata =  (whiteningMatrix * env_concat)';
  fsample = round(1/mode(diff(t)));
  
  save([savefile_hmm 'data'],'hmmdata','whiteningMatrix','fsample','subj_inds')
  save([savefile_hmm 'data'],'subj_inds','-append')
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            I N F E R   H M M                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if todo.infer

    if ~exist('hmmdata','var')
      load([savefile_hmm 'data.mat'])
    end


    hmm = ABhmm_infer(hmmdata,nstates,nreps);
    hmm.whiteningMatrix = whiteningMatrix;
    hmm.fsample = fsample;

    % Save results
    save(savefile_hmm,'hmm','subj_inds')
    HMMresults = [savefile_hmm '.mat'];

    
end; % added MWW

if ~exist('hmm','var') % added MWW  
  load([savefile_hmm '.mat']);
  HMMresults = [savefile_hmm '.mat'];
end

% Compute partial correlation maps
statePath = ABhmm_statepath(hmm);

pcorr = zeros(size(hmm.whiteningMatrix,2),hmm.K);
for subnum = unique(subj_inds)
    
  [~,filestr] = fileparts(BFfiles{subnum});
  disp(['Computing partial correlation maps for ' filestr]);
  load(fullfile(hmmdir,[BFnames{subnum} '_movavg' num2str(windowsize*1000) 'ms']));
  if size(env,3)>1, % More than 1 trial so concatenate them, MWW
    envreshaped=reshape(env,[size(env,1) size(env,2)*size(env,3)]);
  else
    envreshaped=env;
  end;
      
  disp(['ntrials=' num2str(size(env,3))]);
  pcorr = pcorr + ABhmm_regress(statePath(subj_inds==subnum),envreshaped,0,'pcorr');
   
  % get trial-wise state paths for each subject
  
  %subj_inds = [subj_inds,subnum*ones(1,size(env,2))];

  sp_sub=statePath(subj_inds==subnum);

  if size(env,3)>1, % More than 1 trial so concatenate them, MWW
      sub_statePath{subnum}=reshape(sp_sub,[size(env,2),size(env,3)]); % ntpts x ntrials
  end;

end
pcorr = pcorr ./ length(unique(subj_inds));

load(BFfiles{1});
statemaps = fullfile(hmmdir,'/multisub_movavg100ms_pcorr');
statemaps = nii_quicksave(pcorr,statemaps,oat_stage_results.gridstep,2);





