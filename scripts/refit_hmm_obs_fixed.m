%
% Run HMM analysis
%
load([dirs.vars '/prepFiles'], 'prepFiles', 'T');

% Only use first nSubjectsToDo data files
prepFiles = prepFiles(1:nSubjectsToDo);
T = T(1:nSubjectsToDo);
%prepFiles = prepFiles(end-nSubjectsToDo+1:end);
%T = T(end-nSubjectsToDo+1:end);

% HMM options (see setup.m)
options = hmm_options;

% Data preparation (time embedding and PCA) has already been performed
options = rmfield(options, 'embeddedlags');
options = rmfield(options, 'pca');

% Remove stochastic learning options
options = rmfield(options, 'BIGNinitbatch');
options = rmfield(options, 'BIGNbatch');
options = rmfield(options, 'BIGtol');
options = rmfield(options, 'BIGcyc');
options = rmfield(options, 'BIGundertol_tostop');
options = rmfield(options, 'BIGdelay');
options = rmfield(options, 'BIGforgetrate');
options = rmfield(options, 'BIGbase_weights');

% Hold the observation model fixed
options.updateObs = 0;

% Get the fitted hmm
hmm = load([dirs.results '/hmm.mat']);
options.hmm = hmm.hmm;

% Fit HMM
gamma = [];
for i = 1:nSubjectsToDo
    disp(['Fitting subject ' num2str(i)]);
    disp(['===================']);
    [hmm, Gamma, ~, vpath] = hmmmar({prepFiles{i}}, {T{i}}, options);
    gamma = [gamma; Gamma];
end
save([dirs.results '/obs_fixed_gamma.mat'], 'gamma', '-v7.3');

clearvars -except dirs freqRange hmm_options nEmbeddings nStates nSubjectsToDo session;
