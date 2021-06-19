%
% Applies a Hilbert transform to source reconstructed data
%

% Session info
session.name = 'eo'; % eo, vmg, vms, vml
if strcmp(session.name, 'eo')
    session.optPrefix = 'Bffd';
else
    session.optPrefix = 'Reffd';
end

% Directories
dirs.base   = ['/well/woolrich/projects/uk_meg_notts/' session.name];
dirs.srcRec = [dirs.base '/summer21/src_rec'];
dirs.prep   = [dirs.base '/summer21/prepared_data'];

warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir(dirs.prep);

disp('using directories:');
disp(dirs);

%
% Get source reconstructed data
%
fileArray  = dir([dirs.srcRec '/pBFf' session.optPrefix '*.mat']);
nSubjects = length(fileArray);

srcRecFiles = cell(nSubjects, 1);
for i = 1:nSubjects
    srcRecFiles{i} = [fileArray(i).folder '/' fileArray(i).name];
end

%
% Prepare data
%
for i = 1:nSubjects
    fprintf('\nhilbert transforming session %d\n', i);
    fprintf('=================================\n');

    S                     = struct();
    S.D                   = srcRecFiles{i};
    S.winsize             = 1/40; % secs
    S.downsample          = 0;
    S.remove_edge_effects = 1;
    S.prefix              = 'h';

    D = osl_hilbenv(S);
end

% Save as mat files
for i = 1:nSubjects
    disp(['Saving subject ' num2str(i)]);
    [filepath, filename, ext] = fileparts(srcRecFiles{i});
    spmFile = [dirs.srcRec '/' S.prefix filename];
    matFile = [dirs.prep '/subject' num2str(i)];
    read_spm_file(spmFile, matFile);
end
