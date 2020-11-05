%
% Bandpass filter the data and perform beamforming
%
load([dirs.vars '/sessionInfo'], 'nSessions', 'sessionSpmFiles', 'spmLinks');

maskFile  = [osldir '/std_masks/MNI152_T1_8mm_brain.nii.gz'];
mniCoords = osl_mnimask2mnicoords(maskFile);

bfFiles = cell(1, nSessions);
for i = 1:nSessions
    fprintf('\nfiltering and beamforming for session %d\n', i);
    fprintf('==========================================\n');

    % Bandpass filtering
    S      = struct();
    S.band = 'bandpass';
    S.freq = freqRange;
    S.D    = sessionSpmFiles{i};
    D      = spm_eeg_filter(S);

    % Move bandpassed data into beamforming folder
    newD       = D.move([dirs.bf '/']);
    bfFiles{i} = fullfile(newD);
    
    % Beamforming
    S                = struct();
    S.modalities     = {'MEGGRAD'};
    S.timespan       = [0 Inf];
    S.pca_order      = 250;
    S.type           = 'Scalar';
    S.inverse_method = 'beamform';
    S.prefix         = '';
    osl_inverse_model(bfFiles{i}, mniCoords, S);

    % Remove temporary files
    runcmd(['rm -r ' dirs.bf '/*temp*']);
end

save([dirs.vars '/bfFiles'],   'bfFiles');
save([dirs.vars '/maskFile'],  'maskFile');

clear nSessions sessionSpmFiles spmLinks maskFile mniCoords i S D newD bfFiles ans;
