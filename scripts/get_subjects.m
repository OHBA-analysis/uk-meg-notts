%
% Setups up arrays containing paths to subject data
%

% Subjects
subjects  = dir([dirs.rawData '/3*']);
subjects  = sort_nat({subjects.name});
nSubjects = length(subjects);

% Arrays containing paths
ctfFiles = cell(nSubjects, 1);
mriFiles = cell(nSubjects, 1);
posFiles = cell(nSubjects, 1);
for i = 1:nSubjects
    ctfFiles{i} = ...
        add_if_exists([dirs.rawData '/' subjects{i} '/' session.dsFiles]);
    mriFiles{i} = ...
        add_if_exists([dirs.rawData '/' subjects{i} '/' subjects{i} '_CRG.mri']);
    posFiles{i} = add_if_exists([dirs.rawData '/pos_files/' subjects{i} '.pos']);
end

% Remove subjects with missing files
removeSubjectIndices = find(cellfun(@isempty, ctfFiles) | ...
                            cellfun(@isempty, mriFiles) | ...
                            cellfun(@isempty, posFiles));
removeSubjectIndices = flip(removeSubjectIndices); % must reverse because we're deleting
removeSubjectIndices = [removeSubjectIndices; 8];  % remove subject 3018
check_subjects(removeSubjectIndices, nSubjects);
subjects  = remove_subjects(subjects, removeSubjectIndices);
ctfFiles  = remove_subjects(ctfFiles, removeSubjectIndices);
mriFiles  = remove_subjects(mriFiles, removeSubjectIndices);
posFiles  = remove_subjects(posFiles, removeSubjectIndices);
nSubjects = length(subjects);

save([dirs.vars '/subjectInfo'],  'subjects', 'nSubjects');
save([dirs.vars '/rawDataFiles'], 'ctfFiles', 'mriFiles', 'posFiles');

clearvars -except dirs freqRange hmm_options nEmbeddings nStates nSubjectsToDo firstSubject session;

function filePath = add_if_exists(filePath)
    file = dir(filePath);
    if isempty(file)
        filePath = [];
    else
        filePath = [file.folder '/' file.name];
    end
    clear file;
end

function check_subjects(indices, nSubjects)
    nRemove = length(indices);
    if nRemove == nSubjects
        error('all subjects have missing files');
    end
    fprintf('using %d out of %d subjects\n\n', nSubjects - nRemove, nSubjects);
    clear nRemove;
end

function filenameArray = remove_subjects(filenameArray, indices)
    for i = 1:length(indices)
        j = indices(i);
        filenameArray(j) = [];
    end
    clear i j;
end
