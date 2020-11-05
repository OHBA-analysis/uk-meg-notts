%
% Gets the initial covariances used when fitting an HMM
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

% Only train the HMM for 1 cycle
options.cyc = 1;

% Fit HMM
[hmm, Gamma, ~, vpath] = hmmmar(prepFiles, T, options);
save([dirs.results '/hmm_cyc-1.mat'], 'hmm', '-v7.3');

clear prepFiles T options hmm;
