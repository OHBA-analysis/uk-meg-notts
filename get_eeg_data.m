%
% Saves EEG data
%

% Session info
session.name = 'eo'; % eo, vmg, vms, vml
if strcmp(session.name, 'eo')
    session.optPrefix = 'Bffd';
    %session.optPrefix = 'BAffd';
else
    %session.optPrefix = 'Reffd';
    session.optPrefix = 'ffd';
    %session.optPrefix = 'Affd';
end

disp('session info:')
disp(session)

% Directories
dirs.base = ['/well/woolrich/projects/uk_meg_notts/' session.name];
dirs.opt  = [dirs.base '/preproc.opt'];
dirs.eeg  = [dirs.base '/eeg']

disp('using directories:');
disp(dirs);

%
% Get preprocessed data
%
fileArray = dir([dirs.opt '/' session.optPrefix '*.mat']);
nSubjects = length(fileArray);

optSpmFiles = cell(nSubjects, 1);
for i = 1:nSubjects
    optSpmFiles{i} = [fileArray(i).folder '/' fileArray(i).name];
end

%
% Save EEG data
%
for i = 1:nSubjects
    disp(['Saving subject ' num2str(i)]);
    D = spm_eeg_load(optSpmFiles{i});
    inds = [D.indchannel('EEG057'), D.indchannel('EEG058'), D.indchannel('EEG059'), D.indchannel('EEG060')];
    [X, T] = read_spm_file(optSpmFiles{i});
    X = X(:,inds);
    save([dirs.eeg '/subject' num2str(i)], 'X');
end
