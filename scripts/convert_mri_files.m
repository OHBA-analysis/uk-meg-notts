%
% Converts .mri files to .nii files
%
load([dirs.vars '/subjectInfo'],  'subjects', 'nSubjects');
load([dirs.vars '/rawDataFiles'], 'mriFiles');

% MRI conversions need freesurfer
addpath('/home/cgohil/local/freesurfer/matlab')

niiFiles = cell(nSubjects, 1);
for i = 1:nSubjects
    niiFiles{i} = [dirs.nii '/' subjects{i} '/' subjects{i} '_CRG.nii'];
end
save([dirs.vars '/niiFiles'], 'niiFiles');

for i = 1:nSubjects
    fprintf('\nconverting mri file for subject %s\n', subjects{i});
    fprintf('====================================\n');
    runcmd(['mkdir -p ' dirs.nii '/' subjects{i}]);
    convert_mri(mriFiles{i}, niiFiles{i});
end

clear subjects nSubjects mriFiles niiFiles i;
