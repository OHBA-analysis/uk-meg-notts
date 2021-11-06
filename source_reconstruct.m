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
    %session.optPrefix = 'BAffd';
else
    %session.optPrefix = 'Reffd';
    %session.optPrefix = 'ffd';
    session.optPrefix = 'Affd';
end

disp('session info:')
disp(session)

% Data preparation settings
nEmbeddings = 15;

% Directories
dirs.base   = ['/well/woolrich/projects/uk_meg_notts/' session.name];
dirs.opt    = [dirs.base '/preproc.opt'];
%dirs.srcRec = [dirs.base '/natcomms18_ica_' session.optPrefix '/src_rec'];
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
%parcFile  = [osldir '/parcellations/fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'];

filtPrefix = 'f';
bfPrefix   = 'BF';
parcPrefix = 'p';
dipPrefix = 'DSF';

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
    S.freq   = [1 45];
    S.prefix = filtPrefix;
    S.D      = [dirs.srcRec '/' filename ext];

    spm_eeg_filter(S);

    bandpassOptions = S;
    
    % Beamforming
    fprintf('\nBeamforming session %d\n', i);

    S                = struct();
    S.modalities     = {'MEGGRAD'};
    S.timespan       = [0 Inf];
    S.pca_order      = 120;
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
    S.orthogonalisation     = 'symmetric';
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
end

% Parcellated data files
parcFiles = cell(nSubjects, 1);
for i = 1:nSubjects
    [filepath, filename, ext] = fileparts(optSpmFiles{i});
    parcFiles{i} = [dirs.srcRec '/' parcPrefix bfPrefix filtPrefix filename ext];
end


% Fix dipole sign ambiguity
fprintf('\nFixing dipole sign ambiguity\n');

% Add functions needed by this script to MATLAB path
addpath('sign_flipping');

% Establish a good template subject
S = struct();
S.concat = struct();
S.concat.protocol = 'none';
S.concat.embed.do = 1;
S.concat.embed.num_embeddings = nEmbeddings;
S.concat.embed.rectify = false;
S.concat.whiten = 1;
S.concat.normalisation = 'voxelwise';
S.concat.pcadim = -1;
S.netmat_method = @netmat_cov;

state_netmats_cov_preflipped = hmm_full_global_cov(parcFiles, S);

% Assess which subject is the best template
state_netmats = state_netmats_cov_preflipped;

modes       = {'none','abs'};
diag_offset = 15;

metric_global = zeros(length(state_netmats), length(state_netmats), length(modes));
for mm=1:length(modes)
    for subj=1:length(state_netmats)
       for subj2=1:length(state_netmats)
            if subj2 ~= subj
                metric_global(subj, subj2, mm) = matrix_distance_metric( ...
                    state_netmats{subj}.global.netmat_full, ...
                    state_netmats{subj2}.global.netmat_full, ...
                    diag_offset,modes{mm}, ...
                    []);
            end
       end
    end
end

tmp = sum(metric_global(:,:,2), 2);
template_subj = nearest(tmp, median(tmp));

% Perform the sign flip
S = struct();
S.roinets_protocol = parcellationOptions.orthogonalisation;
S.innovations_mar_order = parcellationOptions.innovations_mar_order;
S.Ds = parcFiles;
S.num_iters = 500;
S.prefix = dipPrefix;
S.num_embeddings = nEmbeddings;
S.subj_template = template_subj;
[signflipped_files_out, sign_flip_results] = find_sign_flips(S);

dipoleOptions = S;

% Save source reconstructed data as a mat file
for i = 1:nSubjects
    fprintf('\nSaving session %d\n', i);
    [filepath, filename, ext] = fileparts(optSpmFiles{i});
    srcRecFile  = [dirs.srcRec '/' dipPrefix parcPrefix bfPrefix filtPrefix filename ext];
    matFile = [dirs.srcRec '/subject' num2str(i) '.mat'];
    read_spm_file(srcRecFile, matFile);
end

save([dirs.srcRec '/options'], 'bandpassOptions', 'beamformingOptions', 'parcellationOptions', 'dipoleOptions');
