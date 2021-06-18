%
% Performs source reconstruction on preprocessed data, this includes
% - Beamforming
% - Parcellation
% - Dipole sign fixing
%

% Session info
session.name = 'eo'; % eo, vmg, vms, vml
if strcmp(session.name, 'eo')
    session.optPrefix = 'Bffd';
else
    session.optPrefix = 'Reffd';
end

disp('session info:')
disp(session)

% Directories
dirs.base   = ['/well/woolrich/projects/uk_meg_notts/' session.name];
dirs.opt    = [dirs.base '/preproc.opt'];
dirs.srcRec = [dirs.base '/natcomms18/src_rec'];

warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir(dirs.srcRec);

disp('using directories:');
disp(dirs);

%
% Get preprocessed data
%
fileArray = dir([dirs.opt '/' session.optPrefix '*.mat']);
nSubjects = length(fileArray);

optSpmFiles = cell(nSubjects, 1);
for i = 1:nSubjects
    optSpmFiles{i} = [fileArray(i).folder '/' fileArray(i).name];
end

%
% Source reconstruction
%
maskFile  = [osldir '/std_masks/MNI152_T1_8mm_brain.nii.gz'];
mniCoords = osl_mnimask2mnicoords(maskFile);
parcFile  = [osldir '/parcellations/fmri_d100_parcellation_with_3PCC_ips_reduced_2mm_ss5mm_ds8mm_adj.nii.gz'];

filtPrefix = 'f';
bfPrefix   = 'BF';
parcPrefix = 'p';

matFiles = cell(nSubjects, 1);
T        = cell(nSubjects, 1);
for i = 1:nSubjects

    % Copy preprocessed data into source reconstruction folder
    S   = struct();
    S.D = optSpmFiles{i};
    [filepath, filename, ext] = fileparts(optSpmFiles{i});
    S.outfile = [dirs.srcRec '/' filename ext];

    spm_eeg_copy(S);

    % Bandpass filtering
    fprintf('\nFiltering session %d\n', i);

    S        = struct();
    S.band   = 'bandpass';
    S.freq   = [1 98];
    S.prefix = filtPrefix;
    S.D      = [dirs.srcRec '/' filename ext];

    spm_eeg_filter(S);

    bandpassOptions = S;
    
    % Beamforming
    fprintf('\nBeamforming session %d\n', i);

    S                = struct();
    S.modalities     = {'MEGGRAD'};
    S.timespan       = [0 Inf];
    S.pca_order      = 250;
    S.type           = 'Scalar';
    S.inverse_method = 'beamform';
    S.prefix         = bfPrefix;

    D = [dirs.srcRec '/' filtPrefix filename ext];
    osl_inverse_model(D, mniCoords, S);

    runcmd(['rm -r ' dirs.srcRec '/*temp*']);

    beamformingOptions = S;

    % Parcellation
    fprintf('Parcellating session %d\n', i);

    S                       = struct();
    S.D                     = [dirs.srcRec '/' bfPrefix filtPrefix filename ext];
    S.parcellation          = parcFile;
    S.orthogonalisation     = 'innovations_mar';
    S.innovations_mar_order = 14;
    S.method                = 'spatialBasis';
    S.normalise_voxeldata   = 0;
    S.prefix                = parcPrefix;
    S.maskfname             = maskFile;

    [D, weights, assignments]  = osl_apply_parcellation(S);
    D.parcellation.weights     = weights;
    D.parcellation.assignments = assignments;
    D.save;

    parcellationOptions = S;

    % Save source reconstructed data as a mat file
    fprintf('\nSaving session %d\n', i);

    srcRecFile  = [dirs.srcRec '/' parcPrefix bfPrefix filtPrefix filename ext];
    matFiles{i} = [dirs.srcRec '/subject' num2str(i) '.mat'];
    [~, t]      = read_spm_file(srcRecFile, matFiles{i});
    T{i}        = t;

end

% Fix dipole sign ambiguity
fprintf('\nFixing dipole sign ambiguity\n');

S         = struct();
S.maxlag  = 7;
s.verbose = 1;

flips = findflip(matFiles, T, S);
flipdata(matFiles, T, flips, [], 1);

dipoleOptions = S;

save([dirs.srcRec '/options'], 'bandpassOptions', 'beamformingOptions', 'parcellationOptions', 'dipoleOptions');

clear;
