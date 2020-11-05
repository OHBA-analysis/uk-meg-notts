%
% Manually check and correct coregistration
%
load_opt;

for i = 1:nSubjects
    D = spm_eeg_load(opt.results.spm_files{i});
    if S.use_rhino
        rhino_manual(D);
    else
        spm_eeg_inv_checkdatareg(D);
        %spm_eeg_inv_checkmeshes(D, 1);
    end
    waitfor(gcf)
end

% Re-run forward model if position edited in rhino_manual
for i = 1:nSubjects
    S.D = opt.results.spm_files{i};
    S.forward_meg = 'MEG Local Spheres';
    osl_forward_model(S);
end

clear opt i D S;
