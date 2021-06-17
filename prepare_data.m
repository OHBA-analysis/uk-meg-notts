%
% Data preparation for an HMM, this includes
% - Time embedding
% - PCA
% - Standardisation
%

% Session info
session.name = 'eo'; % eo, vmg, vms, vml

% Directories
dirs.base   = ['/well/woolrich/projects/uk_meg_notts/' session.name];
dirs.srcRec = [dirs.base '/natcomms18/src_rec'];
dirs.prep   = [dirs.base '/natcomms18/prepared_data'];

warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir(dirs.prep);

disp('using directories:');
disp(dirs);

% HMM options
options.embeddedlags   = -7:7;
options.pca            = 80;
options.standardise    = 1;
options.standardise_pc = 1;

% These options are needed for validation but aren't used
options.K         = 2;
options.BIGNbatch = 1;

%
% Get source reconstructed data
%
fileArray  = dir([dirs.srcRec '/subject*.mat']);
nSubjects = length(fileArray);

srcRecFiles = cell(nSubjects, 1);
T           = cell(nSubjects, 1);
for i = 1:nSubjects
    srcRecFiles{i} = [fileArray(i).folder '/' fileArray(i).name];
    data = load(srcRecFiles{i});
    T{i} = data.T;
end

%
% Prepare data
%

% Account for discontinuities in the time series
%T_all = T; 

% Pretend the entire time series is continuous
T_all = cell(nSubjects, 1);
for i = 1:nSubjects
    T_all{i} = sum(T{i});
end

% Check options, this adds defaults to options
[options, srcRecFiles] = checkoptions(options, srcRecFiles, T_all, 0);

%
% Data preparation from hmmmar.m
%
disp('Calculating PCA');

% get PCA pre-embedded loadings
if length(options.pca_spatial) > 1 || (options.pca_spatial > 0 ...
   && options.pca_spatial ~= 1)
    if ~isfield(options,'As')
        options.As = highdim_pca(srcRecFiles,T_all,options.pca_spatial,...
            0,options.standardise,...
            options.onpower,0,options.detrend,...
            options.filter,options.leakagecorr,options.Fs);
    end
    options.pca_spatial = size(options.As,2);
else
    options.As = [];
end

if length(options.embeddedlags) > 1
    elmsg = '(embedded)';
else
    elmsg = '';
end

% get PCA loadings
if length(options.pca) > 1 || (options.pca > 0 && options.pca ~= 1) || ...
        isfield(options,'A')
    if ~isfield(options,'A')
        [options.A,~,e] = highdim_pca(srcRecFiles,T_all,options.pca,...
            options.embeddedlags,options.standardise,...
            options.onpower,options.varimax,options.detrend,...
            options.filter,options.leakagecorr,options.Fs,options.As);
        options.pca = size(options.A,2);
        if options.verbose
            if options.varimax
                fprintf('Working in PCA/Varimax %s space, with %d components.\n',...
                        elmsg,options.pca)
                fprintf('(explained variance = %1f)  \n',e(options.pca))
            else
                fprintf('Working in PCA %s space, with %d components.\n',...
                        elmsg,options.pca)
                fprintf('(explained variance = %1f)\n',...
                        e(options.pca))
            end
        end
    end
    options.ndim = size(options.A,2);
    options.S = ones(options.ndim);
    options.Sind = formindexes(options.orders,options.S);
    if ~options.zeromean
        options.Sind = [true(1,size(options.Sind,2)); options.Sind];
    end
else
    options.As = [];
end

if isfield(options,'A') && ~isempty(options.A)
    options.ndim = size(options.A,2);
elseif isfield(options,'As') && ~isempty(options.As)
    options.ndim = size(options.As,2);
else
    X = loadfile(srcRecFiles{1},T_all{1},options);
    options.ndim = size(X,2);
end

if options.pcamar > 0 && ~isfield(options,'B')
    % PCA on the predictors of the MAR regression, per lag:
    % X_t = \sum_i X_t-i * B_i * W_i + e
    options.B = pcamar_decomp(srcRecFiles,T_all,options);
end

if options.pcapred > 0 && ~isfield(options,'V')
    % PCA on the predictors of the MAR regression, together: 
    % Y = X * V * W + e, where X contains all the lagged predictors
    % So, unlike B, V draws from the temporal dimension and not only spatial
    options.V = pcapred_decomp(srcRecFiles,T_all,options);
end

% Save prepared data for each subject to file
for i = 1:nSubjects
    disp(['Preparing data for subject ' num2str(i)]);
    data = cell(1, 1);
    data{1} = srcRecFiles{i};
    T = cell(1, 1);
    T{1} = T_all{i};
    [X,~,~,T] = loadfile(data, T, options); % loadfile does the time embedding and PCA
    prepFile = [dirs.prep '/subject' num2str(i) '.mat'];
    save(prepFile, 'X', 'T');
end

% Save the PCA weights
W = options.A;
save([dirs.prep '/pca_components'], 'W');

clear;
