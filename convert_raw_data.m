%
% Converts raw data files to compatible formats for OSL
%

% Session info
session.name = 'vml'; % eo, vmg, vms, vml
switch session.name
    case 'eo'
        session.dsFiles = '*Eyes_Open*.ds';
    case 'vmg'
        session.dsFiles = '*VisMotor_Gamma*.ds';
    case 'vms'
        session.dsFiles = '*VisMotor_Short*.ds';
    case 'vml'
        session.dsFiles = '*VisMotor_Long*.ds';
end

disp('session info:')
disp(session)

% Directories
dirs.data = '/well/woolrich/projects/uk_meg_notts/raw_data';
dirs.raw  = ['/well/woolrich/projects/uk_meg_notts/' session.name '/raw'];

warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir(dirs.raw);

disp('using directories:');
disp(dirs);

%
% Get subjects
%
subjects  = dir([dirs.data '/3*']);
subjects  = sort_nat({subjects.name});
nSubjects = length(subjects);

% Arrays containing paths
ctfFiles = cell(nSubjects, 1);
mriFiles = cell(nSubjects, 1);
posFiles = cell(nSubjects, 1);
for i = 1:nSubjects
    ctfFiles{i} = add_if_exists([dirs.data '/' subjects{i} '/' session.dsFiles]);
    mriFiles{i} = add_if_exists([dirs.data '/' subjects{i} '/' subjects{i} '_CRG.mri']);
    posFiles{i} = add_if_exists([dirs.data '/pos_files/' subjects{i} '.pos']);
end

% Remove subjects with missing files
removeSubjectIndices = find(cellfun(@isempty, ctfFiles) | cellfun(@isempty, mriFiles) | cellfun(@isempty, posFiles));
removeSubjectIndices = flip(removeSubjectIndices); % must reverse because we're deleting
removeSubjectIndices = [removeSubjectIndices; 8];  % remove subject 3018
check_subjects(removeSubjectIndices, nSubjects);
subjects  = remove_subjects(subjects, removeSubjectIndices);
ctfFiles  = remove_subjects(ctfFiles, removeSubjectIndices);
mriFiles  = remove_subjects(mriFiles, removeSubjectIndices);
posFiles  = remove_subjects(posFiles, removeSubjectIndices);
nSubjects = length(subjects);

% Create directories for each subject's raw data
subjectDirs = cell(nSubjects, 1);
for i = 1:nSubjects
    subjectDirs{i} = [dirs.raw '/' subjects{i}];
    mkdir(subjectDirs{i});
end

%
% Convert ds files to SPM files
%
for i = 1:nSubjects
    fprintf('\nconverting ds file for subject %s\n', subjects{i});
    fprintf('===================================\n');

    S = struct();
    S.outfile = [subjectDirs{i} '/' subjects{i} '.mat'];
    D = osl_import(ctfFiles{i}, S);

    % Fiducials
    posFile = fopen(posFiles{i});
    fiducialData = textscan(posFile, '%s %f %f %f');
    fclose(posFile);

    fiducialIndices = [find(strcmpi(fiducialData{1}, 'nasion')) find(strcmpi(fiducialData{1}, 'left')) find(strcmpi(fiducialData{1}, 'right'))];

    newFiducials.fid.pnt   = [fiducialData{2}(fiducialIndices) fiducialData{3}(fiducialIndices) fiducialData{4}(fiducialIndices)] * 10;
    newFiducials.fid.label = {'nas'; 'lpa'; 'rpa'};

    % Headshape
    headShapeIndices  = setdiff(2:length(fiducialData{1}), fiducialIndices);

    newFiducials.pnt  = [fiducialData{2}(headShapeIndices) fiducialData{3}(headShapeIndices) fiducialData{4}(headShapeIndices)] * 10;
    newFiducials.unit = 'mm';

    D = fiducials(D, newFiducials);
    D.save;
end

%
% Convert structural MRI files
%
addpath('/well/woolrich/projects/software/freesurfer/matlab')

for i = 1:nSubjects
    fprintf('\nconverting mri file for subject %s\n', subjects{i});
    fprintf('====================================\n');

    niiFile = [subjectDirs{i} '/' subjects{i} '.nii'];
    convert_mri(mriFiles{i}, niiFile);
end

clear;

%
% Helper functions
%
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

function convert_mri(mriFileName, niftiFileName)
    [fileDir, fileName, fileExt] = fileparts(niftiFileName);
    
    % Begin conversion
    mri = ft_read_mri(mriFileName);
    mri.anatomy = flip(mri.anatomy,1);
    mri.anatomy = flip(mri.anatomy,2);
    mri.anatomy = flip(mri.anatomy,3);
    mri.anatomy = flip(mri.anatomy,1);

    % write out
    ft_write_mri(niftiFileName, mri.anatomy, 'dataformat','nifti');

    % This bit added by MWW to make sure hdr info is correct
    nii = load_untouch_nii(niftiFileName);

    [pth nm ext]=fileparts(niftiFileName);

    % We assume qform and sform are the same, and correspond to the
    % transform from voxel coords to native scanner space.
    % For some reason rhino does not cope well if native scanner
    % space coords are too different from MNI. So to help with this
    % we change the offset to get the native scanner space coords
    % looking a bit more like MNI (note that the chosen native 
    % scanner space coordinate system is arbitrary wrt what rhino 
    % needs to do). (presumably it is a local minima
    % issue in rhino - and ideally should be fixed in there).
    % We need to make sure the qformcode and sformcode is 1.
    % As the output from edit_header
    % sets an invalid code which rhino will not like.
    % See convert_mri for more on this
    offset=[128 -128 -90];
    [in,out] = edit_header(niftiFileName, fullfile(pth,nm), nii.hdr.dime.pixdim(2:4), offset, -1);

    % add correct sform code.
    % the output from edit_header has an invalid qform and
    % seemingly correct sform (we can check sform is working as
    % when it is viewed in fsleyes, all of the labels seem to be
    % correct (although can not be sure about left-right or the offset) ).
    % However, the sform code is set to 4. For rhino to be happy
    % the sform code needs to be set to 1. So we manually enforce
    % this now. (Note that this call also seems to set the qform to be
    % the same as the sform):
    % fslorient -setsformcode 1 imagename
    runcmd(['fslorient -setsformcode 1 ' niftiFileName]);
end

function [in,out] = edit_header(input_file, output_file, pixdim, offset, qfac)
    if nargin < 5, qfac=-1; end
    if nargin < 4, offset=zeros(1,3); end
    if nargin < 3, pixdim=ones(1,3); end

    assert( ischar(input_file) && ischar(output_file), 'First and second inputs must filenames.' );
    assert( isnumeric(pixdim) && numel(pixdim)==3 && all(pixdim>0), 'Third input should be a 1x3 positive vector of pixdims.' );
    assert( isnumeric(offset) && numel(offset)==3, 'Fourth input should be a 1x3 translation vector.' );
    assert( isnumeric(qfac) && isscalar(qfac), 'Fifth input should be a scalar q-factor.' );

    % load input nifti file
    in = load_untouch_nii(input_file);

    % information needed for solving the orientation
    dat = in.hdr.hist;
    dim = in.hdr.dime;

    % for affine transformations only, the simplest is to use the s-form
    dat.qform_code = 0;
    dat.sform_code = 4;
    %dat.magic = 'n+1';

    % set pixdims
    dim.pixdim(1) = qfac;
    dim.pixdim([2 3 4]) = pixdim;

    % set s-form
    dat.srow_x = [ qfac*pixdim(1)         0         0 offset(1) ];
    dat.srow_y = [              0 pixdim(2)         0 offset(2) ];
    dat.srow_z = [              0         0 pixdim(3) offset(3) ];

    % set corresponding q-form too just to be sure
    dat.qoffset_x = offset(1);
    dat.qoffset_y = offset(2);
    dat.qoffset_z = offset(3);

    dat.quatern_b = 0;
    dat.quatern_c = 0;
    dat.quatern_d = 0;

    % save output file with modified headers
    out = in;
    out.hdr.hist = dat;
    out.hdr.dime = dim;

    fprintf('Saving output file: "%s"\n',output_file);
    save_untouch_nii( out, output_file );
end
