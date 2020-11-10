%
% Runs the OSL proprocessing tool (OPT)
%
load([dirs.vars '/subjectInfo'], 'nSubjects');
load([dirs.vars '/niiFiles'],    'niiFiles')
load([dirs.vars '/spmFiles'],    'spmFiles');

optName = [session.name '240520'];
save([dirs.vars '/optName'], 'optName');

% Label artefact channels
for i = 1:nSubjects
    D = spm_eeg_load(spmFiles{i});
    D = D.chantype(find(strcmp(D.chanlabels, 'EEG057')), 'EOG1');
    D = D.chantype(find(strcmp(D.chanlabels, 'EEG058')), 'EOG2');
    D = D.chantype(find(strcmp(D.chanlabels, 'EEG059')), 'ECG');
    D = D.chantype(find(strcmp(D.chanlabels, 'EEG060')), 'EMG');
    D.save;
end

% Settings
opt           = struct();
opt.dirname   = [dirs.data '/' optName];
opt.datatype  = 'ctf';
opt.spm_files = spmFiles;

opt.downsample.do   = 1;
opt.downsample.freq = 250;

% Highpass and notch filtering
opt.highpass.do     = 1;
opt.highpass.cutoff = 0.1;
opt.mains.do        = 0;

% AFRICA denoising
%opt.africa.do          = 1;
opt.africa.todo.ica    = 0;
opt.africa.todo.ident  = 0;
opt.africa.todo.remove = 0;
opt.africa.ident.func  = @identify_artefactual_components_auto;
opt.africa.ident.kurtosis_wthresh       = 0.2;
opt.africa.ident.max_num_artefact_comps = 2;
opt.africa.precompute_topos             = 1;

nEventTypes = length(session.eventTypes);
opt.bad_segments.do = ~(nEventTypes > 0);
opt.bad_segments.event_significance = 0.05;

% Coregistration
opt.coreg.do           = true;
opt.coreg.mri          = niiFiles;
opt.coreg.use_rhino    = 1;
opt.coreg.useheadshape = 1;
opt.coreg.forward_meg  = 'MEG Local Spheres';

% Epoch settings
opt.epoch.do = (nEventTypes > 0);
if opt.epoch.do
    for i = 1:nEventTypes
        opt.epoch.trialdef(i).conditionlabel = session.eventTypes{i};
        opt.epoch.trialdef(i).eventtype      = session.eventTypes{i};
        opt.epoch.trialdef(i).eventvalue     = 1;
        opt.epoch.time_range                 = session.eventTimeRange;
    end
    opt.outliers.do                  = true;
    opt.outliers.outlier_measure_fns = {'std'};
    %opt.outliers.event_significance = 0.05;
else
    opt.outliers.do = 0;
end

% Run OPT
%opt = opt_consolidate_results(opt);
opt = osl_run_opt(opt);

clear nSubjects niiFiles spmFiles optName i D S nEventTypes opt ans;
