%
% Computes spectra at the group level and per subject
%
load([dirs.vars '/preprocFiles'], 'preprocFiles', 'T');
load([dirs.results '/hmm'], 'hmm')

preprocFiles = preprocFiles(1:nSubjectsToDo);
T = T(1:nSubjectsToDo);
%preprocFiles = preprocFiles(end-nSubjectsToDo+1:end);
%T = T(end-nSubjectsToDo+1:end);

optionsMt = struct();
optionsMt.Fs = 250; % Sampling rate
optionsMt.fpass = freqRange; % band of frequency you're interested in
optionsMt.tapers = [4 7]; % taper specification
optionsMt.p = 0; % interval of confidence  
optionsMt.win = 500; % multitaper window
optionsMt.to_do = [1 0]; % turn off pdc
optionsMt.order = 0;
optionsMt.embeddedlags = -round(nEmbeddings/2):round(nEmbeddings/2);

% average
disp(['Computing spectra for all subjects']);
fitMt = hmmspectramt(preprocFiles, T, hmm.Gamma, optionsMt);

% per subject
fitMtSubject = cell(length(preprocFiles), 1);
d = length(optionsMt.embeddedlags) - 1;
acc = 0;
for n = 1:length(preprocFiles)
    disp(['Computing spectra for subject ' num2str(n)]);
    load(preprocFiles{n}, 'X', 'T');
    gamma = hmm.Gamma(acc + (1:(sum(T)-length(T)*d)),:);
    acc = acc + size(gamma, 1);
    fitMtSubject{n} = hmmspectramt(X, T, gamma, optionsMt);
    fitMtSubject{n}.state = rmfield(fitMtSubject{n}.state, 'ipsd');
    fitMtSubject{n}.state = rmfield(fitMtSubject{n}.state, 'pcoh');
    fitMtSubject{n}.state = rmfield(fitMtSubject{n}.state, 'phase');
end
save([dirs.results '/fitMt'], 'fitMtSubject', 'fitMt');

clear preprocFiles X T hmm optionsMt fitMt d acc n gamma fitMtSubject;
