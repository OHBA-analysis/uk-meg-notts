%
% Prepares data for training an HMM
%
load([dirs.vars '/preprocFiles'], 'preprocFiles', 'T');
T_all = T;
nSubjects = length(T_all);

% HMM options (see setup.m)
options = hmm_options;
[options,preprocFiles] = checkoptions(options,preprocFiles,T_all,0);

%--------------------------------%
% Data preparation from hmmmar.m %
%--------------------------------%
disp('Calculating PCA');

% get PCA pre-embedded loadings
if length(options.pca_spatial) > 1 || (options.pca_spatial > 0 ...
   && options.pca_spatial ~= 1)
    if ~isfield(options,'As')
        options.As = highdim_pca(preprocFiles,T_all,options.pca_spatial,...
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
        [options.A,~,e] = highdim_pca(preprocFiles,T_all,options.pca,...
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
    X = loadfile(preprocFiles{1},T_all{1},options);
    options.ndim = size(X,2);
end

if options.pcamar > 0 && ~isfield(options,'B')
    % PCA on the predictors of the MAR regression, per lag:
    % X_t = \sum_i X_t-i * B_i * W_i + e
    options.B = pcamar_decomp(preprocFiles,T_all,options);
end

if options.pcapred > 0 && ~isfield(options,'V')
    % PCA on the predictors of the MAR regression, together: 
    % Y = X * V * W + e, where X contains all the lagged predictors
    % So, unlike B, V draws from the temporal dimension and not only spatial
    options.V = pcapred_decomp(preprocFiles,T_all,options);
end

% Save prepared data for each subject to file
prepFiles = cell(nSubjects, 1);
T_all_new = cell(nSubjects, 1);
for i = 1:nSubjects
    disp(['Preparing data for subject ' num2str(i)]);
    data = cell(1, 1);
    data{1} = preprocFiles{i};
    T = cell(1, 1);
    %T{1} = T_all{i}; % account for discontinuities in the time series
    T{1} = sum(T_all{i}); % pretend the entire time series is continuous
    [X,~,~,T] = loadfile(data, T, options); % loadfile does the time embedding and PCA
    prepFiles{i} = [dirs.prepData '/subject' num2str(i) '.mat'];
    save(prepFiles{i}, 'X', 'T');
    T_all_new{i} = T;
end
T = T_all_new;
save([dirs.vars '/prepFiles'], 'prepFiles', 'T', '-v7.3');

% Save the PCA weights
W = options.A;
save([dirs.prepData '/pca_weights'], 'W', '-v7.3');

clear preprocFiles X T T_all options ans data e elmsg i nSubjects W;
