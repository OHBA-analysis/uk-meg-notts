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

% Fit HMM
[hmm, Gamma, ~, vpath] = hmmmar(prepFiles, T, options);
hmm.Gamma = Gamma;
hmm.vpath = vpath;
hmm.T = T;
save([dirs.results '/hmm.mat'], 'hmm', '-v7.3');

clearvars -except dirs freqRange hmm_options nEmbeddings nStates nSubjectsToDo session;
