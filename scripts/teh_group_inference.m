%
% Run HMM analysis using Mark's method,
% which calls osl-core/osl_hmm_toolbox/teh_groupinference_parcels.m
%
load([dirs.vars '/signFlippedFiles'], 'signFlippedFiles');
load([dirs.vars '/parcSettings'], 'parcS');

% Struct to hold HMM options
hmmoptions = struct();

% Dimensionality reduction
D = spm_eeg_load(signFlippedFiles{1});
nparcels = size(D.parcellation.weights, 2);
hmmoptions.prepare.pcadim = 80;

% Time Embedding 
hmmoptions.prepare.embed.do = 1;
hmmoptions.prepare.embed.num_embeddings = nEmbeddings;

% Number of states
hmmoptions.hmm.K = 12;

hmmoptions.prepare.normalisation = 'voxelwise';
hmmoptions.prepare.whiten        = 1;
hmmoptions.prepare.savePCmaps    = 0;
%hmmoptions.prepare.max_ntpts     = 40000;

hmmoptions.hmm.dynamic_model_type = 'hmm';
%hmmoptions.hmm.dynamic_model_type = 'vbrt';

% To get initial covariances
%hmmoptions.hmm.cyc = 1

hmmoptions.hmm.initcyc   = 60;
hmmoptions.hmm.initrep   = 5;
hmmoptions.hmm.big       = 1;
hmmoptions.hmm.BIGNbatch = 5;

hmmoptions.output.method             = 'pcorr';
hmmoptions.output.use_parcel_weights = 0;
hmmoptions.output.assignment         = 'hard';

% Setup filenames
hmmoptions.prepare.filename = ['hmm_' parcS.prefix ...
                               'pcdim-' num2str(hmmoptions.prepare.pcadim) ...
                               '_' num2str(hmmoptions.prepare.normalisation) ...
                               '_embed-' num2str(hmmoptions.prepare.embed.num_embeddings)];

hmmoptions.hmm.filename     = [hmmoptions.prepare.filename ...
                               '_K-' num2str(hmmoptions.hmm.K)...
                               '_big-' num2str(hmmoptions.hmm.big)...
                               '_dyn_model-' hmmoptions.hmm.dynamic_model_type];

hmmoptions.output.filename  = [hmmoptions.hmm.filename '_output'];

% Run HMM
hmmoptions.todo.prepare  = 1;
hmmoptions.todo.hmm      = 1;
hmmoptions.todo.output   = 1;                                                        

hmmoptions.hmmdir = dirs.hmm;

teh_groupinference_parcels(signFlippedFiles, hmmoptions);

clear D ans hmmoptions i nSessions nparcels parcS signFlippedFiles ...
      signflipped_files_out tmp;
