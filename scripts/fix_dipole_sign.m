%
% Fixes dipole sign ambiguity
%
load([dirs.vars '/sessionInfo'],  'nSessions');
load([dirs.vars '/parcFiles'],    'parcFiles');
load([dirs.vars '/parcSettings'], 'parcS');

% Add functions needed by this script to MATLAB path
addpath([dirs.work '/scripts/sign_flipping']);

% Establish a good template subject
S = struct();
S.concat = struct();
S.concat.protocol = 'none';
S.concat.embed.do = 1;
S.concat.embed.num_embeddings = nEmbeddings;
S.concat.embed.rectify = false;
S.concat.whiten = 1;
S.concat.normalisation = 'voxelwise';
S.concat.pcadim = -1;
S.netmat_method = @netmat_cov;

state_netmats_cov_preflipped = hmm_full_global_cov(parcFiles, S);

% Assess which subject is the best template
state_netmats = state_netmats_cov_preflipped;

modes       = {'none','abs'};
diag_offset = 15;

metric_global = zeros(length(state_netmats), length(state_netmats), length(modes));
for mm=1:length(modes)
    for subj=1:length(state_netmats)
       for subj2=1:length(state_netmats)
            if subj2 ~= subj
                metric_global(subj, subj2, mm) = matrix_distance_metric( ...
                    state_netmats{subj}.global.netmat_full, ...
                    state_netmats{subj2}.global.netmat_full, ...
                    diag_offset,modes{mm}, ...
                    []);
            end
       end
    end
end

tmp = sum(metric_global(:,:,2), 2);
template_subj = nearest(tmp, median(tmp));

% Perform the sign flip
S = struct();
S.roinets_protocol = parcS.orthogonalisation;
S.innovations_mar_order = parcS.innovations_mar_order;
S.Ds = parcFiles;
S.num_iters = 500;
S.prefix = 'sfold_';
S.num_embeddings = nEmbeddings;
S.subj_template = template_subj;
[signflipped_files_out, sign_flip_results] = find_sign_flips(S);

sign_flip_results.signflipped_files = signflipped_files_out;
sign_flip_results.energies_group = mean(sign_flip_results.energies, 2);
sign_flip_results.energies = sign_flip_results.energies(1:20:end,:);
sign_flip_results.subj_template_no = template_subj;
sign_flip_results.freq_range = freqRange;

% Paths to files containing sign flipped data
signFlippedFiles = cell(nSessions, 1);
for i = 1:nSessions
    signFlippedFiles{i} = prefix(parcFiles{i}, S.prefix);
end
save([dirs.vars '/signFlippedFiles'], 'signFlippedFiles');

clear parcFiles S state_netmats_cov_preflipped state_netmats modes ...
      diag_offset metric_global mm subj subj2 template_subj sign_flip_results ...
      sign_flipped_files_out signFlippedFiles tmp parcS;
