%
% Fits an HMM to prepared data
%

% Session info
session.name = 'vml'; % eo, vmg, vms, vml

% Directories
dirs.base    = ['/well/woolrich/projects/uk_meg_notts/' session.name];
dirs.prep    = [dirs.base '/natcomms18/prepared_data'];
dirs.results = [dirs.base '/natcomms18/results/Subj1-1_K-6'];

warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir(dirs.results);

disp('using directories:');
disp(dirs);

% Get the first and last subject and number of states from the folder name
splitDir     = split(dirs.results, '/');
params       = regexp(splitDir{end}, '\d*', 'Match');
firstSubject = str2num(params{1});
lastSubject  = str2num(params{2});
nStates      = str2num(params{3});
nSubjects    = lastSubject - firstSubject + 1;


% Fit to time-embedded, PCA data
options          = struct();
options.K        = nStates;
options.order    = 0;
options.zeromean = 1;
options.covtype  = 'full';
options.Fs       = 250;
options.onpower  = 0;
options.inittype = 'HMM-MAR';
options.cyc      = 100;
options.initcyc  = 10;
options.initrep  = 3;
options.verbose  = 1;

% Fit amplitude envelope data
%options          = struct();
%options.K        = nStates;
%options.order    = 0;
%options.covtype  = 'full';
%options.zeromean = 0;
%options.onpower  = 0;
%options.Fs       = 250;
%options.verbose  = 1;

% Stochastic learning options
options.BIGNinitbatch      = 4; % must be less than the number of subjects
options.BIGNbatch          = 4; % must be less than the number of subjects
options.BIGtol             = 1e-7;
options.BIGcyc             = 500;
options.BIGundertol_tostop = 5;
options.BIGdelay           = 5;
options.BIGforgetrate      = 0.7;
options.BIGbase_weights    = 0.9;

options.useParallel = false;

%
% Get prepared data
%
prepFiles = cell(nSubjects, 1);
prepT     = cell(nSubjects, 1);
for i = 1:nSubjects
    prepFiles{i} = [dirs.prep '/subject' num2str(i) '.mat'];
    data = load(prepFiles{i});
    prepT{i} = data.T;
end

%
% Fit an HMM
%
[hmm, Gamma, ~, vpath] = hmmmar(prepFiles, prepT, options);
hmm.Gamma   = Gamma;
hmm.vpath   = vpath;
hmm.T       = prepT;
hmm.options = options;

save([dirs.results '/hmm.mat'], 'hmm', '-v7.3');
