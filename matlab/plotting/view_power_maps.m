mapFile = '/Users/gohil/data/uk_meg_notts/eo/natcomms18/results/Subj1-55_K-8/wb_power_comp1.nii.gz';

workbenchDir = '/Users/gohil/Applications/workbench/bin_macosx64';
path1 = getenv('PATH');
if ~contains(path1, workbenchDir)
    path1 = [path1 ':' workbenchDir];
    setenv('PATH', path1);
end

osl_render4D(mapFile, 'interptype', 'trilinear', 'visualise', true);