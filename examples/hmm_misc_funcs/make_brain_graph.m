%make_brain_graph
% plots the network matrix

figureStem = 'controls brainnet';

% get network matrix
% resultsDir = '/Users/gilesc/data/Henry/OIL-output/analysis-for-parcellation-paper/resting-EO/used_in_paper/ROInets-eo-fmriBasis100reduced-betterReg-closest';
resultsDir = pwd;
networkFile = fullfile(resultsDir, 'ROInetworks_correlation_mats.mat'); %_adj.mat');
load(networkFile);
nROIs = ROInets.cols(correlationMats{1}.correlation);

% get node co-ordinates
parcellationDirectory = '/Users/gilesc/data/papers/parcellation-orthogonalisation/parcellations/fMRI-parcellations/groupPCA_d1000_d100.dr';
% parcelFile = fullfile(parcellationDirectory, 'fmri_d100_parcellation_with_PCC_reduced_2mm_ds8mm.nii.gz');
parcelFile = fullfile(parcellationDirectory, 'fmri_d100_parcellation_with_PCC_tighterMay15_2mm_ds6mm.nii.gz');
spatialRes = 6;
spatialMap = nii_quickread(parcelFile, spatialRes);
mni_coords = find_ROI_centres(spatialMap, spatialRes, 0);

BAND_NAMES = {'Alpha', 'Beta'};
iBand = 2;

% threshold network matrix
graph = correlationMats{iBand}.groupEnvPartialCorrelationRegularized_z;
graphThresh = correlationMats{iBand}.falseDiscoveryRate.zThreshold.groupEnvPartialCorrelationRegularized_z;
graph = -X;
nROIs = rows(X);
graphThresh = 0.05;
ignoreInds = abs(graph) < 4.2; %graphThresh;
graph(ignoreInds) = NaN;

% plot
colorLims = [-8 8];
sphereCols = repmat([30 144 255]/255, nROIs, 1);
edgeLims   = [4.2 8];

figure('Color', 'w', 'Name', [figureStem ' ' BAND_NAMES{iBand} ' band']);
osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], 0, sphereCols, edgeLims);
% 
% % visualise path lengths
% irrelevantInds = tril(true(nROIs), 1);
% graph(irrelevantInds) = NaN;
% [i, j] = find(~isnan(graph));
% for p = length(i):-1:1,
%     pathLength(p) = sqrt(sum((mni_coords(i(p), :) - mni_coords(j(p), :)).^2));
% end
% 
% figure('Color', 'w', 'Name', [figureStem ' path lengths ' BAND_NAMES{iBand} ' band']);
% nBins = 50;
% GC_histogram(pathLength, nBins);
% xlim([0 120]);