%
% Calculate a surface mesh for each subject using the preprocessed data
%

% Session info
session.name = 'eo'; % eo, vmg, vms, vml
if strcmp(session.name, 'eo')
    session.optPrefix = 'Bffd';
else
    %session.optPrefix = 'Reffd';
    session.optPrefix = 'ffd';
end

disp('session info:')
disp(session)

% Directories
dirs.base   = ['/well/woolrich/projects/uk_meg_notts/' session.name];
dirs.opt    = [dirs.base '/preproc.opt'];
dirs.surfMesh = [dirs.base '/australia21/surf_mesh'];
dirs.cwd = pwd;

warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir(dirs.surfMesh);

disp('using directories:');
disp(dirs);

% Copy preprocessed data
fileArray = dir([dirs.opt '/' session.optPrefix '*.mat']);
nSubjects = length(fileArray);

%
% Calculate surface mesh
%
for i = 1:nSubjects

    % Copy and load SPM file containing the preprocessed data
    copyfile([fileArray(i).folder '/' fileArray(i).name],[dirs.surfMesh '/' fileArray(i).name]);
    copyfile([fileArray(i).folder '/' fileArray(i).name(1:end-4) '.dat'],[dirs.surfMesh '/' fileArray(i).name(1:end-4) '.dat']);

    spmFile = [dirs.surfMesh '/' fileArray(i).name];
    D = spm_eeg_load(spmFile);

    S = struct();
    [~,tn] = fileparts(tempname(D.path));
    S.dirname = fullfile(D.path,sprintf('osl_bf_temp_%s',tn(3:3+7)));
    if ~exist(S.dirname,'dir')
        mkdir(S.dirname);
    end

    matlabbatch{1}.spm.tools.beamforming.data.dir = {S.dirname};
    matlabbatch{1}.spm.tools.beamforming.data.D = {D.fullfile};
    matlabbatch{1}.spm.tools.beamforming.data.val = 1;
    matlabbatch{1}.spm.tools.beamforming.data.gradsource = 'inv';
    matlabbatch{1}.spm.tools.beamforming.data.space = 'MNI-aligned';
    matlabbatch{1}.spm.tools.beamforming.data.overwrite = 0;
    matlabbatch{2}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{2}.spm.tools.beamforming.sources.reduce_rank = [2 3];
    matlabbatch{2}.spm.tools.beamforming.sources.keep3d = 1;
    matlabbatch{2}.spm.tools.beamforming.sources.plugin.mesh.orient = 'original';
    matlabbatch{2}.spm.tools.beamforming.sources.plugin.mesh.fdownsample = 1;
    matlabbatch{2}.spm.tools.beamforming.sources.plugin.mesh.flip = false;
    matlabbatch{2}.spm.tools.beamforming.sources.visualise = 0;

    matlabbatch{3}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{3}.spm.tools.beamforming.features.whatconditions.all = 1;
    if D.ntrials == 1
        matlabbatch{3}.spm.tools.beamforming.features.plugin.contcov = struct([]);
    else
        matlabbatch{3}.spm.tools.beamforming.features.plugin.cov = struct([]);
    end

    matlabbatch{3}.spm.tools.beamforming.features.woi = [-Inf Inf]; % needs to be in msecs for bf_features
    matlabbatch{3}.spm.tools.beamforming.features.modality = {'MEGGRAD'};
    matlabbatch{3}.spm.tools.beamforming.features.fuse = 'no';
    matlabbatch{3}.spm.tools.beamforming.features.regularisation.manual.lambda = 0;
    matlabbatch{3}.spm.tools.beamforming.features.bootstrap = false;

    matlabbatch{4}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv_multicov.pca_order = D.nchannels;
    matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv_multicov.type = 'scalar';
    matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv_multicov.bilateral = 0;

    matlabbatch{5}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{5}.spm.tools.beamforming.output.plugin.montage_osl.normalise = 'both';

    matlabbatch{6}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg_osl.prefix = '';

    spm_jobman('run',matlabbatch)

    % Restore the working directory and delete temporary directory
    cd(dirs.cwd);
    delete([S.dirname '/*']);
    rmdir(S.dirname); 
end
