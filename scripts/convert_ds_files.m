%
% Converts .ds files to SPM MEEG objects
%
load([dirs.vars '/subjectInfo'],  'subjects', 'nSubjects');
load([dirs.vars '/rawDataFiles'], 'ctfFiles', 'posFiles');

spmFiles = cell(nSubjects, 1);
for i = 1:nSubjects
    spmFiles{i} = [dirs.spm '/' subjects{i} '.mat'];
end
save([dirs.vars '/spmFiles'], 'spmFiles');

for i = 1:nSubjects
    fprintf('\nconverting ds file for subject %s\n', subjects{i});
    fprintf('===================================\n');

    S = struct();
    S.outfile = spmFiles{i};
    D = osl_import(ctfFiles{i}, S);

    % Fiducials
    posFile = fopen(posFiles{i});
    fiducialData = textscan(posFile, '%s %f %f %f');
    fclose(posFile);

    fiducialIndices = [find(strcmpi(fiducialData{1}, 'nasion')) ...
                       find(strcmpi(fiducialData{1}, 'left')) ...
                       find(strcmpi(fiducialData{1}, 'right'))];

    newFiducials.fid.pnt   = [fiducialData{2}(fiducialIndices) ...
                              fiducialData{3}(fiducialIndices) ...
                              fiducialData{4}(fiducialIndices)] * 10;
    newFiducials.fid.label = {'nas'; 'lpa'; 'rpa'};
    
    % Headshape
    headShapeIndices  = setdiff(2:length(fiducialData{1}), fiducialIndices);

    newFiducials.pnt  = [fiducialData{2}(headShapeIndices) ...
                         fiducialData{3}(headShapeIndices) ...
                         fiducialData{4}(headShapeIndices)] * 10;
    newFiducials.unit = 'mm';

    D = fiducials(D, newFiducials);
    D.save;
end

clear subjects nSubjects ctfFiles posFiles spmFiles i S D posFile ...
      fiducialData fiducialIndices headShapeIndices newFiducials ans;
