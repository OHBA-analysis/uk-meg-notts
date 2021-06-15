%
% Performs a Hilbert transformation to get amplitude envelops
%
load([dirs.vars '/sessionInfo'], 'nSessions');
load([dirs.vars '/parcFiles'], 'parcFiles');

% Hilbert transform
for i = 1:nSessions
    fprintf('\nhilbert transforming session %d\n', i);
    fprintf('=================================\n');

    S                     = struct();
    S.D                   = parcFiles{i};
    S.winsize             = 1/40; % secs
    S.downsample          = 0;
    S.remove_edge_effects = 1;
    S.prefix              = 'h';

    D = osl_hilbenv(S);
end

% Save as standard matlab files
matDir = [dirs.bf '/envelope_mat_files'];
mkdir(matDir)
for i = 1:nSessions
    disp(['Saving subject ' num2str(i)]);
    [filepath, spmFilename, ext] = fileparts(parcFiles{i});
    spmFile = [dirs.bf '/' S.prefix spmFilename];
    matFile = [matDir '/subject' num2str(i)];
    read_spm_file(spmFile, matFile);
end

clearvars -except dirs freqRange hmm_options nEmbeddings nStates nSubjectsToDo session;
