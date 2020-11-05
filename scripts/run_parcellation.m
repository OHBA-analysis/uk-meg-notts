%
% Performs parcellation and orthogonalisation
%
% Leakage correction can be done using orthogonalisation:
% = 'symmetric': the method described in Colclough et al (2015)
% = 'innovations_mar' :the method described in Pasqual-Marqui et al (2017)
%    we also need to specify the option 'innovations_mar_order'
% = 'none'
%
load([dirs.vars '/sessionInfo'], 'nSessions');
load([dirs.vars '/bfFiles'],     'bfFiles');
load([dirs.vars '/maskFile'],    'maskFile');

% Diego's settings:
%parcName          = '42ROI_separatePCC';
%orthogonalisation = 'innovations_mar';

%parcFile = [osldir '/parcellations' ...
%           '/fmri_d100_parcellation_with_3PCC_ips_reduced_2mm_ss5mm_ds8mm_adj.nii.gz'];
%parcPrefix = ['42ROI_separatePCC_' ...
%              orthogonalisation num2str(innovations_mar_order) '_'];

% Mark's settings
orthogonalisation = 'symmetric';
parcFile   = [osldir '/parcellations' ...
              '/fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'];
parcPrefix = ['giles_' orthogonalisation '_'];

% Paths to files containing parcellated data
parcFiles = cell(nSessions, 1);
for i = 1:nSessions
    parcFiles{i} = prefix(bfFiles{i}, parcPrefix);
end

for i = 1:nSessions
    fprintf('\nparcellation for session %d\n', i);
    fprintf('=============================\n');

    parcS                       = struct();
    parcS.D                     = bfFiles{i};
    parcS.parcellation          = parcFile;
    parcS.orthogonalisation     = orthogonalisation;
    parcS.innovations_mar_order = 14;
    parcS.method                = 'spatialBasis';
    parcS.normalise_voxeldata   = 0;
    parcS.prefix                = parcPrefix;
    parcS.maskfname             = maskFile;
    [parcD, parcWeights, parcAssignments] = osl_apply_parcellation(parcS);
    parcD.parcellation.weights     = parcWeights;
    parcD.parcellation.assignments = parcAssignments;
    parcD.save; % save SPM file
end

save([dirs.vars '/parcSettings'], 'parcS');
save([dirs.vars '/parcFiles'],    'parcFiles');

clear nSessions bfFiles maskFile parcName orthogonalisation i parcS parcFile ...
      parcPrefix parcFiles parcD parcWeights parcAssignments ans;
