%
% Splits sessions into separate SPM files
%
load_opt;
load([dirs.vars '/optName'], 'optName');

nSessions       = length(opt.results.spm_files);
sessionSpmFiles = cell(1, nSessions);
spmLinks        = cell(nSessions, 2);
for i = 1:nSessions
    sessionSpmFiles{i} = [dirs.sess '/' optName '_session' num2str(i)];

    S         = struct();
    S.D       = opt.results.spm_files{i};
    S.outfile = sessionSpmFiles{i};
    spm_eeg_copy(S);

    % Log the file each copy corresponds to
    spmLinks{i,1} = S.D;
    spmLinks{i,2} = S.outfile;
end
save([dirs.vars '/sessionInfo'], 'nSessions', 'sessionSpmFiles', 'spmLinks');

clear opt nSessions sessionSpmFiles spmLinks i optName S ans;
