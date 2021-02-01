%
% Settings and directories
%
addpath([pwd '/scripts']);

% Session specific settings
session.name = 'eo'; % eo, vmg, vms, vml
switch session.name
    case 'eo'
        session.dsFiles        = '*Eyes_Open*.ds';
        session.spmFile        = '_eo.mat';
        session.eventTypes     = [];
        session.eventTimeRange = [-10 10];
    case 'vmg'
        session.dsFiles        = '*VisMotor_Gamma*.ds';
        session.spmFile        = '_vmg.mat';
        session.eventTypes     = {'Abduction', 'Stim_Onset'};
        session.eventTimeRange = [-4 4];
    case 'vms'
        session.dsFiles        = '*VisMotor_Short*.ds';
        session.spmFile        = '_vms.mat';
        session.eventTypes     = {'Abduction', 'Stim_Onset'};
        session.eventTimeRange = [-4 4];
    case 'vml'
        session.dsFiles        = '*VisMotor_Long*.ds';
        session.spmFile        = '_vml.mat';
        session.eventTypes     = {'Abduction', 'Stim_Onset'};
        session.eventTimeRange = [-4 4];
end

% HMM settings
freqRange   = [1 45];
nEmbeddings = 15; % must be an odd number
nStates     = 12;
nSubjectsToDo = 45;

hmm_options = struct();
hmm_options.order = 0;
hmm_options.zeromean = 1;
hmm_options.covtype = 'full';
hmm_options.embeddedlags = -floor(nEmbeddings/2):floor(nEmbeddings/2);
hmm_options.pca = 80;
hmm_options.K = nStates;
hmm_options.Fs = 250;
hmm_options.verbose = 1;
hmm_options.onpower = 0;
hmm_options.standardise = 1;
hmm_options.standardise_pc = hmm_options.standardise;
hmm_options.inittype = 'HMM-MAR';
hmm_options.cyc = 100;
hmm_options.initcyc = 10;
hmm_options.initrep = 3;

% Stochastic learning options
hmm_options.BIGNinitbatch = 5; % must be less than the number of subjects
hmm_options.BIGNbatch = 5; % must be less than the number of subjects
hmm_options.BIGtol = 1e-7;
hmm_options.BIGcyc = 500;
hmm_options.BIGundertol_tostop = 5;
hmm_options.BIGdelay = 5;
hmm_options.BIGforgetrate = 0.7;
hmm_options.BIGbase_weights = 0.9;

hmm_options.useParallel = false;

% Directories
dirs.work    = pwd;
dirs.vars    = [dirs.work '/variables'];

dirs.rawData     = '/well/woolrich/shared/uk_meg_notts/raw_data';
dirs.data        = ['/well/woolrich/shared/uk_meg_notts/' session.name];
dirs.spm         = [dirs.data '/spm'];
dirs.sess        = [dirs.data '/sessions'];
dirs.nii         = [dirs.data '/nii'];
dirs.bf          = [dirs.data '/bf'];
dirs.preprocData = [dirs.data '/preproc_data'];
dirs.prepData    = [dirs.data '/prepared_data'];
dirs.results     = [dirs.data '/results/nSubjects-' num2str(nSubjectsToDo) '_K-' num2str(hmm_options.K)];
%dirs.results     = [dirs.data '/results/nSubjects-' num2str(nSubjectsToDo) '_K-' num2str(hmm_options.K) '_hold-out'];

warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir(dirs.spm)
mkdir(dirs.sess)
mkdir(dirs.nii)
mkdir(dirs.bf)
mkdir(dirs.preprocData)
mkdir(dirs.prepData)
mkdir(dirs.results)

disp('using directories:');
disp(dirs);
