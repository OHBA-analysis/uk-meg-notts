load('/Users/gohil/data/uk_meg_notts/eo/natcomms18/results/Subj1-55_K-8/wb_coh_comp1');

K = size(coh,1);
nROIs = size(coh,2);
M = zeros(nROIs); 
for k = 1:K
    M = M + squeeze(abs(coh(k,:,:))) / K;
end
coh = squeeze(coh(7,:,:)) - M;

tmp=triu(coh);
inds2=find(tmp>1e-10);
data=tmp(inds2);

if length(data)>0
    S2=[];
    S2.data=squash(data);
    S2.do_fischer_xform=false;
    S2.pvalue_th = 0.01/length(S2.data);
    S2.do_plots=1;
    graph_ggm=teh_graph_gmm_fit(S2);
    th = graph_ggm.orig_th;
    graph=coh;
else
    graph = nan(nparcels);
    inds2 = find(ones(nparcels));
end
graph(graph<th)=NaN;

colorLims = [th-1 th*1.2]; 
sphereCols = repmat([30 144 255]/255, nROIs, 1); 
edgeLims = [4 8];

% specify parcellation:
%parcellation_file = [osldir '/parcellations/fmri_d100_parcellation_with_3PCC_ips_reduced_2mm_ss5mm_ds8mm_adj.nii.gz'];
parcellation_file = [osldir '/parcellations/fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'];
spatialRes = 8; 
spatialMap = nii.quickread(parcellation_file, spatialRes); 
mni_coords = find_ROI_centres(spatialMap, spatialRes, 0); 

figure('Color', 'w','Position',[547 100 577 453]);
ax(1) = axes('Position',[0 0.5 0.5 0.5]);
ax(2) = axes('Position',[0.55 0.5 0.5 0.5]);
ax(3) = axes('Position',[0.27 0.1 0.5 0.5]);
viewZ = {[270,0],[-270,0],[0,90]};
for iplot=1:3
    axes(ax(iplot))
    osl_braingraph(graph, colorLims, repmat(0.5,nROIs,1), [0 1], mni_coords, [], 0, sphereCols, edgeLims); 
    view(ax(iplot),viewZ{iplot})
    colorbar('hide') 
end