%
% Gets maps in NIFTI and CIFTI format
%
load([dirs.vars '/maskFile'], 'maskFile');
load([dirs.vars '/parcSettings'], 'parcS');
load([dirs.results '/fitMt'], 'fitMtGroupFactWb', 'fitMtGroupFact4b');
load([dirs.results '/hmm'], 'hmm');

% Setup workbench
workbenchDir = '/well/woolrich/projects/software/workbench/bin_rh_linux64';
path1 = getenv('PATH');
if ~contains(path1, workbenchDir)
    path1 = [path1 ':' workbenchDir];
    setenv('PATH', path1);
end

[mask, res, xform] = nii.load(maskFile);
parcFile = parcS.parcellation;
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
osl_render4D([mapFile '.nii.gz'], 'savedir', mapFile, ...
             'interptype', 'trilinear', 'visualise', false)

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
    osl_render4D([mapFile '.nii.gz'], 'savedir', mapFile, ...
                 'interptype', 'trilinear', 'visualise', false)
end

clearvars -except dirs freqRange hmm_options nEmbeddings nStates nSubjectsToDo session;
