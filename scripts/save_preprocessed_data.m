%
% Save data to MATLAB files for the HMM fitting
%
load([dirs.vars '/signFlippedFiles'], 'signFlippedFiles');

nSubjects = length(signFlippedFiles);

% Prepare data files for the HMM
preprocFiles = cell(nSubjects, 1);
T = cell(nSubjects, 1);
for i = 1:nSubjects
    disp(['Saving preprocessed data for subject ' num2str(i)]);
    preprocFiles{i} = [dirs.preprocData '/subject' num2str(i) '.mat'];
    [~, t] = read_spm_file(signFlippedFiles{i}, preprocFiles{i});
    T{i} = t;
end

save([dirs.vars '/preprocFiles'], 'preprocFiles', 'T');

clear signFlippedFiles nSubjects preprocFiles T i t;
