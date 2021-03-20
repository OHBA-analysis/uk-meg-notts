%
% Saves data to a matlab file and fixes dipole sign ambiguity
%
load([dirs.vars '/parcFiles'], 'parcFiles');

nSubjects = length(parcFiles);

% Save as matlab files
preprocFiles = cell(nSubjects, 1);
T = cell(nSubjects, 1);
for i = 1:nSubjects
    disp(['Saving preprocessed data for subject ' num2str(i)]);
    preprocFiles{i} = [dirs.preprocData '/subject' num2str(i) '.mat'];
    [~, t] = read_spm_file(parcFiles{i}, preprocFiles{i});
    T{i} = t;
end

save([dirs.vars '/preprocFiles'], 'preprocFiles', 'T');

% Fix dipole sign ambguity
disp('Fixing dipole sign ambiguity');

options = struct();
options.maxlag = floor(nEmbeddings / 2);
options.verbose = 1;

flips = findflip(preprocFiles, T, options);
flipdata(preprocFiles, T, flips);

clearvars -except dirs freqRange hmm_options nEmbeddings nStates nSubjectsToDo session;
