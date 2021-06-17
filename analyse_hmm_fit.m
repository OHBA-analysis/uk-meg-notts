%
% Analyses an HMM fit
%

% Session info
session.name = 'eo'; % eo, vmg, vms, vml

% Directories
dirs.base    = ['/well/woolrich/projects/uk_meg_notts/' session.name];
dirs.srcRec  = [dirs.base '/natcomms18/src_rec'];
dirs.results = [dirs.base '/natcomms18/results/Subj1-55_K-12'];

disp('using directories:');
disp(dirs);

%
% Get the HMM fit
%
hmm = load([dirs.results '/hmm.mat'], 'hmm');
hmm = hmm.hmm;

% Get the first and last subject and number of states from the folder name
splitDir     = split(dirs.results, '/');
params       = regexp(splitDir{end}, '\d*', 'Match');
firstSubject = str2num(params{1});
lastSubject  = str2num(params{2});
nStates      = str2num(params{3});

%
% Get source reconstructed data
%
fileArray = dir([dirs.srcRec '/subject*.mat']);
nSubjects = length(fileArray);

srcRecFiles = cell(nSubjects, 1);
srcRecT     = cell(nSubjects, 1);
for i = 1:nSubjects
    srcRecFiles{i} = [fileArray(i).folder '/' fileArray(i).name];
    data = load(srcRecFiles{i});
    srcRecT{i} = data.T;
end

srcRecFiles = srcRecFiles(firstSubject:lastSubject);
srcRecT     = srcRecT(firstSubject:lastSubject);

%
% Compute spectra
%
options              = struct();
options.Fs           = 250;      % Sampling rate
options.fpass        = [1 45];   % band of frequency you're interested in
options.tapers       = [4 7];    % taper specification
options.p            = 0;        % interval of confidence  
options.win          = 500;      % multitaper window
options.to_do        = [1 0];    % turn off pdc
options.order        = 0;
options.embeddedlags = -7:7;

% Group level
fprintf('Computing spectra for all subjects\n');

fitMt = hmmspectramt(srcRecFiles, srcRecT, hmm.Gamma, options);

% Subject level
fitMtSubject = cell(nSubjects, 1);
d   = length(options.embeddedlags) - 1;
acc = 0;
for i = 1:nSubjects
    disp(['Computing spectra for subject ' num2str(i)]);
    load(srcRecFiles{i}, 'X', 'T');
    gamma = hmm.Gamma(acc + (1:(sum(T)-length(T)*d)),:);
    acc = acc + size(gamma, 1);
    fitMtSubject{i} = hmmspectramt(X, T, gamma, options);
    fitMtSubject{i}.state = rmfield(fitMtSubject{i}.state, 'ipsd');
    fitMtSubject{i}.state = rmfield(fitMtSubject{i}.state, 'pcoh');
    fitMtSubject{i}.state = rmfield(fitMtSubject{i}.state, 'phase');
end
save([dirs.results '/fitMt'], 'fitMtSubject', 'fitMt');

%
% Spectral factorisation
%
fprintf('\nComputing spectral factorisation\n');

% Get the three bands depicted in the 2018 Nature Comms paper,
% the 4th is essentially capturing noise
options       = struct();
options.Ncomp = 4;
options.Base  = 'coh';

[fitMtGroupFact4b, spProfiles4b, fitMtSubjectFact4b] = spectdecompose(fitMtSubject, options);

save([dirs.results '/fitMt'], 'fitMtSubjectFact4b', 'fitMtGroupFact4b', 'spProfiles4b', '-append')

% Get the wideband maps (the second is capturing noise)
options.Ncomp = 2;

[fitMtGroupFactWb, spProfileWb, fitMtSubjectFactWb] = spectdecompose(fitMtSubject, options);

save([dirs.results '/fitMt'], 'fitMtSubjectFactWb', 'fitMtGroupFactWb', 'spProfileWb', '-append')

% Check if the spectral profiles make sense, if not you may want to repeat
fig = figure;
subplot(1,2,1);
plot(spProfiles4b, 'LineWidth', 2);
subplot(1,2,2);
plot(spProfileWb, 'LineWidth', 2);
saveas(fig, [dirs.results '/profiles.png']);

%
% Calculate power maps
%
fprintf('\nCalculating power maps\n');

% Setup workbench
workbenchDir = '/well/woolrich/projects/software/workbench/bin_rh_linux64';
path1 = getenv('PATH');
if ~contains(path1, workbenchDir)
    path1 = [path1 ':' workbenchDir];
    setenv('PATH', path1);
end

% Mask and parcellation file
maskFile = [osldir '/std_masks/MNI152_T1_8mm_brain.nii.gz'];
parcFile = [osldir '/parcellations/fmri_d100_parcellation_with_3PCC_ips_reduced_2mm_ss5mm_ds8mm_adj.nii.gz'];

[mask, res, xform] = nii.load(maskFile);
spatialMap = parcellation(parcFile);
spatialMap = spatialMap.to_matrix(spatialMap.weight_mask);

% Normalise the parcels to have comparable weights 
for j = 1:size(spatialMap, 2)
    spatialMap(:,j) =  spatialMap(:,j) / max(spatialMap(:,j));
end

% Wideband
mapFile = [dirs.results '/state_maps_wideband'];
map = zeros(size(spatialMap,1), length(fitMtGroupFactWb.state));
for k = 1:hmm.K
    psd = diag(squeeze(fitMtGroupFactWb.state(k).psd(1,:,:)));
    map(:,k) = spatialMap * psd;
end
% Center by voxel across states
map = map - repmat(mean(map,2), 1, size(map,2));
mapFile = [dirs.results '/state_maps_wideband'];
nii.save(matrix2vols(map, mask), res, xform, mapFile);
osl_render4D([mapFile '.nii.gz'], 'savedir', mapFile, 'interptype', 'trilinear', 'visualise', false)

% Per frequency band
for fr = 1:3
    mapFile = [dirs.results '/state_maps_band' num2str(fr)];
    map = zeros(size(spatialMap,1), length(fitMtGroupFact4b.state));
    for k = 1:hmm.K
        psd = diag(squeeze(fitMtGroupFact4b.state(k).psd(fr,:,:)));
        map(:,k) = spatialMap * psd;
    end
    % Center by voxel across states
    map = map - repmat(mean(map,2), 1, size(map,2));
    mapFile = [dirs.results '/state_maps_band' num2str(fr)];
    nii.save(matrix2vols(map, mask), res, xform, mapFile);
    osl_render4D([mapFile '.nii.gz'], 'savedir', mapFile, 'interptype', 'trilinear', 'visualise', false)
end

clear;
