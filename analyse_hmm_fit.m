%
% Analyses an HMM fit
%

% Session info
session.name = 'eo'; % eo, vmg, vms, vml

% Directories
dirs.base    = ['/well/woolrich/projects/uk_meg_notts/' session.name];
dirs.srcRec  = [dirs.base '/natcomms18/src_rec'];
dirs.results = [dirs.base '/natcomms18/results/Subj1-55_K-8'];

disp('using directories:');
disp(dirs);

%
% Get the HMM fit
%
hmm = load([dirs.results '/hmm.mat'], 'hmm');
hmm = hmm.hmm;

% Get the first and last subject and number of states from the folder name
splitDir     = split(dirs.results, '/');
params       = regexp(splitDir{end}, '\d*', 'Match');
firstSubject = str2num(params{1});
lastSubject  = str2num(params{2});
nStates      = str2num(params{3});
nSubjects    = lastSubject - firstSubject + 1;

%
% Get source reconstructed data
%
srcRecFiles = cell(nSubjects, 1);
srcRecT     = cell(nSubjects, 1);
for i = 1:nSubjects
    srcRecFiles{i} = [dirs.srcRec '/subject' num2str(i) '.mat'];
    data = load(srcRecFiles{i});
    srcRecT{i} = data.T;
end

%
% Compute spectra
%
options              = struct();
options.Fs           = 250;      % Sampling rate
options.fpass        = [1 45];   % band of frequency you're interested in
options.tapers       = [4 7];    % taper specification
options.p            = 0;        % interval of confidence  
options.win          = 500;      % multitaper window
options.to_do        = [1 0];    % turn off pdc
options.order        = 0;
options.embeddedlags = -7:7;

% Group level
fprintf('Computing spectra for all subjects\n');

fitMt = hmmspectramt(srcRecFiles, srcRecT, hmm.Gamma, options);
save([dirs.results '/fitMt'], 'fitMt', '-v7.3');

% Subject level
fitMtSubject = cell(nSubjects, 1);
d   = length(options.embeddedlags) - 1;
acc = 0;
for i = 1:nSubjects
    disp(['Computing spectra for subject ' num2str(i)]);
    load(srcRecFiles{i}, 'X', 'T');
    gamma = hmm.Gamma(acc + (1:(sum(T)-length(T)*d)),:);
    acc = acc + size(gamma, 1);
    fitMtSubject{i} = hmmspectramt(X, T, gamma, options);
    fitMtSubject{i}.state = rmfield(fitMtSubject{i}.state, 'ipsd');
    fitMtSubject{i}.state = rmfield(fitMtSubject{i}.state, 'pcoh');
    fitMtSubject{i}.state = rmfield(fitMtSubject{i}.state, 'phase');
end
save([dirs.results '/fitMt'], 'fitMtSubject', '-append');

%
% Spectral factorisation
%
fprintf('\nComputing spectral factorisation\n');

% Get the three bands depicted in the 2018 Nature Comms paper,
% the 4th is essentially capturing noise
options       = struct();
options.Ncomp = 4;
options.Base  = 'coh';
options.Niterations = 1;

[fitMtGroupFact4b, spProfiles4b, fitMtSubjectFact4b] = spectdecompose(fitMtSubject, options);

save([dirs.results '/fitMtFact'], '-v7.3', 'fitMtSubjectFact4b', 'fitMtGroupFact4b', 'spProfiles4b');

% Get the wideband maps (the second is capturing noise)
options.Ncomp = 2;

[fitMtGroupFactWb, spProfileWb, fitMtSubjectFactWb] = spectdecompose(fitMtSubject, options);

save([dirs.results '/fitMtFact'], 'fitMtSubjectFactWb', 'fitMtGroupFactWb', 'spProfileWb', '-append');

% Check if the spectral profiles make sense, if not you may want to repeat
fig = figure;
subplot(1,2,1);
plot(spProfiles4b, 'LineWidth', 2);
subplot(1,2,2);
plot(spProfileWb, 'LineWidth', 2);
saveas(fig, [dirs.results '/profiles.png']);
