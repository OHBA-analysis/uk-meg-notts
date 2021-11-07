%
% Save power and coherence maps
%

% Session info
session.name = 'eo'; % eo, vmg, vms, vml

% Directories
dirs.base    = ['/well/woolrich/projects/uk_meg_notts/' session.name];
dirs.results = [dirs.base '/natcomms18/results/Subj1-55_K-8'];

disp('using directories:');
disp(dirs);

%
% Load spectral decomposition
%
load([dirs.results '/fitMtFact'], 'fitMtGroupFactWb', 'fitMtGroupFact4b');

nWbComp = size(fitMtGroupFactWb.state(1).psd, 1);
n4bComp = size(fitMtGroupFact4b.state(1).psd, 1);
nStates = length(fitMtGroupFactWb.state);
nNodes = size(fitMtGroupFactWb.state(1).psd, 2);

%
% Save maps
%

% Source reconstruction files
parcFile  = [osldir '/parcellations/fmri_d100_parcellation_with_3PCC_ips_reduced_2mm_ss5mm_ds8mm_adj.nii.gz'];
%parcFile  = [osldir '/parcellations/fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'];

maskFile = [osldir '/std_masks/MNI152_T1_8mm_brain.nii.gz'];
[mask, res, xform] = nii.load(maskFile);

% Loop through wideband spectral components
for i=1:nWbComp-1

    % Power maps
    spatialMap = parcellation(parcFile);
    spatialMap = spatialMap.to_matrix(spatialMap.weight_mask);
    for j=1:size(spatialMap,2)
        spatialMap(:,j) =  spatialMap(:,j) / max(spatialMap(:,j));
    end
    map = zeros(size(spatialMap,1),nStates);
    for k = 1:nStates
        psd = diag(squeeze(fitMtGroupFactWb.state(k).psd(i,:,:)));
        map(:,k) = spatialMap * psd;
    end
    map = map - repmat(mean(map,2), 1, size(map,2));
    nii.save(matrix2vols(map,mask), res, xform, [dirs.results '/wb_power_comp' num2str(i)]);
    
    % Coherence maps
    coh = zeros(nStates,nNodes,nNodes);
    for k = 1:nStates
        coh(k,:,:) = squeeze(fitMtGroupFactWb.state(k).coh(i,:,:));
    end
    save([dirs.results '/wb_coh_comp' num2str(i)], 'coh');

end

% Loop through narrowband spectral components
for i=1:n4bComp-1

    % Power maps
    spatialMap = parcellation(parcFile);
    spatialMap = spatialMap.to_matrix(spatialMap.weight_mask);
    for j=1:size(spatialMap,2)
        spatialMap(:,j) =  spatialMap(:,j) / max(spatialMap(:,j));
    end
    map = zeros(size(spatialMap,1),nStates);
    for k = 1:nStates
        psd = diag(squeeze(fitMtGroupFact4b.state(k).psd(i,:,:)));
        map(:,k) = spatialMap * psd;
    end
    map = map - repmat(mean(map,2), 1, size(map,2));
    nii.save(matrix2vols(map,mask), res, xform, [dirs.results '/nb_power_comp' num2str(i)]);
    
    % Coherence maps
    coh = zeros(nStates,nNodes,nNodes);
    for k = 1:nStates
        coh(k,:,:) = squeeze(fitMtGroupFact4b.state(k).coh(i,:,:));
    end
    save([dirs.results '/nb_coh_comp' num2str(i)], 'coh');

end

