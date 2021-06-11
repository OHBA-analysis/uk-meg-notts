%
% Run HMM analysis
%
load([dirs.vars '/prepFiles'], 'prepFiles', 'T');

prepFiles = prepFiles(firstSubject:firstSubject + nSubjectsToDo - 1);
T = T(firstSubject:firstSubject + nSubjectsToDo - 1);

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

clearvars -except dirs freqRange hmm_options nEmbeddings nStates nSubjectsToDo firstSubject session;
