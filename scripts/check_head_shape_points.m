%
% Manually check head shape points
%
load([dirs.spm '/spmFiles'], 'spmFiles');

for i = 1:nSubjects
    D = spm_eeg_load(spmFiles{i});
    osl_edit_fid(D)
    waitfor(gcf)
end

clearvars -except dirs freqRange hmm_options nEmbeddings nStates nSubjectsToDo firstSubject session;
