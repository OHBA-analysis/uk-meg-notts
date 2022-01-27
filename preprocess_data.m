%
% Preprocesses raw MEG data using OPT
%

% Session info
session.name = 'eo'; % eo, vmg, vms, vml
switch session.name
    case 'eo'
        session.eventTypes     = [];
    case 'vmg'
        session.eventTypes     = {'Abduction', 'Stim_Onset'};
        session.eventTimeRange = [-10 10];
    case 'vms'
        session.eventTypes     = {'Abduction', 'Stim_Onset'};
        session.eventTimeRange = [-4 4];
    case 'vml'
        session.eventTypes     = {'Abduction', 'Stim_Onset'};
        session.eventTimeRange = [-10 10];
end

disp('session info:')
disp(session)

% Directories
dirs.base    = ['/well/woolrich/projects/uk_meg_notts/' session.name];
dirs.raw     = [dirs.base '/raw'];
dirs.preproc = [dirs.base '/preproc'];

disp('using directories:');
disp(dirs);

%
% Get files needed for OPT
%
subjects  = dir(dirs.raw);
subjects  = subjects(~ismember({subjects.name},{'.','..'}));
nSubjects = length(subjects);

spmFiles = cell(nSubjects, 1);
niiFiles = cell(nSubjects, 1);
for i = 1:nSubjects
    spmFiles{i} = [subjects(i).folder '/' subjects(i).name '/' subjects(i).name '.mat'];
    niiFiles{i} = [subjects(i).folder '/' subjects(i).name '/' subjects(i).name '.nii'];
end

%
% OPT
%

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
opt.dirname   = dirs.preproc;
opt.datatype  = 'ctf';
opt.spm_files = spmFiles;

opt.downsample.do   = 1;
opt.downsample.freq = 250;

% Highpass and notch filtering
opt.highpass.do     = 1;
opt.highpass.cutoff = 0.1;
opt.mains.do        = 1;

% AFRICA denoising
opt.africa.do          = 0;
opt.africa.todo.ica    = 0;
opt.africa.todo.ident  = 'auto';
opt.africa.todo.remove = 1;
opt.africa.ident.func  = @identify_artefactual_components_auto;
opt.africa.ident.kurtosis_wthresh       = 0.2;
%opt.africa.ident.max_num_artefact_comps = 2;
opt.africa.precompute_topos             = 1;
opt.africa.ident.artefact_chans         = {'EOG1', 'EOG2', 'ECG', 'EMG'};

% Bad segment removal
nEventTypes = length(session.eventTypes);
opt.bad_segments.do = ~(nEventTypes > 0);
%opt.bad_segments.event_significance = 0.05;

% Coregistration
opt.coreg.do           = 1;
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
opt = osl_run_opt(opt);
