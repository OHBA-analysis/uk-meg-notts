%
% Performs automatic spectral factorisation to find spectrally defined networks
% Manual inspection of profiles (profiles.png) is advised.
%
load([dirs.vars '/preprocFiles'], 'preprocFiles');
load([dirs.results '/fitMt'], 'fitMtSubject');
load([dirs.results '/hmm'], 'hmm');

preprocFiles = preprocFiles(1:nSubjectsToDo);
%preprocFiles = preprocFiles(end-nSubjectsToDo+1:end);

% Get the three bands depicted in the paper (the 4th is essentially capturing noise)
optionsFact = struct();
optionsFact.Ncomp = 4;
optionsFact.Base = 'coh';
[fitMtGroupFact4b, spProfiles4b, fitMtSubjectFact4b] = ...
    spectdecompose(fitMtSubject, optionsFact);
save([dirs.results '/fitMt'], ...
     'fitMtSubjectFact4b', 'fitMtGroupFact4b', 'spProfiles4b', '-append')

% Get the wideband maps (the second is capturing noise)
optionsFact.Ncomp = 2;
[fitMtGroupFactWb, spProfileWb, fitMtSubjectFactWb] = ...
    spectdecompose(fitMtSubject, optionsFact);
save([dirs.results '/fitMt'], ...
     'fitMtSubjectFactWb', 'fitMtGroupFactWb', 'spProfileWb', '-append')

% Check if the spectral profiles make sense, if not you may want to repeat
fig = figure;
subplot(1,2,1);
plot(spProfiles4b, 'LineWidth', 2);
subplot(1,2,2);
plot(spProfileWb, 'LineWidth', 2);
saveas(fig, [dirs.results '/profiles.png']);

% Do statistical testing on the spectral information
%fitMtSubjectFact1d = cell(length(preprocFiles), 1);
%for n = 1:length(preprocFiles)
%    disp(['Performing statistical testing for subject ' num2str(n)]);
%    fitMtSubjectFact1d{n} = struct();
%    fitMtSubjectFact1d{n}.state = struct();
%    for k = 1:hmm.K % we don't care about the second component
%        fitMtSubjectFact1d{n}.state(k).psd = fitMtSubjectFactWb{n}.state(k).psd(1,:,:);
%        fitMtSubjectFact1d{n}.state(k).coh = fitMtSubjectFactWb{n}.state(k).coh(1,:,:);
%    end
%end
%testsSpectra = specttest(fitMtSubjectFact1d, 5000, 1, 1);
%significantSpectra = spectsignificance(testsSpectra, 0.01);
%save([dirs.results '/fitMt'], 'testsSpectra', 'significantSpectra', '-append')

clearvars -except dirs freqRange hmm_options nEmbeddings nStates nSubjectsToDo session;
