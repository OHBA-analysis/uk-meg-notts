% PLI for HCP

leakageCorrected = true;

%% Directories
dataDir    = '~/scratch/data/hcp-MEG/';
sourceDir  = fullfile(dataDir, 'source-timecourses/uncorrected/');
anatomyDir = fullfile(dataDir, 'anatomy/');
resultsDir = fullfile(dataDir, 'network-analysis/PLI-uncorrected/allCorrected');
parcDir    = fullfile(dataDir, 'parcellations/');
metaDir    = fullfile(dataDir, 'metadata');
ROInets.make_directory(resultsDir);

%% Get file list
deltaDataFiles = osl_filelist(sourceDir, '*1-4Hz_ROI_timecourses.mat');
thetaDataFiles = osl_filelist(sourceDir, '*4-8Hz_ROI_timecourses.mat');
alphaDataFiles = osl_filelist(sourceDir, '*8-13Hz_ROI_timecourses.mat');
betaDataFiles  = osl_filelist(sourceDir, '*13-30Hz_ROI_timecourses.mat');
gammaDataFiles = osl_filelist(sourceDir, '*30-48Hz_ROI_timecourses.mat');

nFiles    = length(deltaDataFiles);


%% Pipeline
%parpool(3)
windowLength = 2; % s
Fs           = 300; % Hz

for iFile = length(deltaDataFiles):-1:1,
    % get file info
    [~,deltaFiles{iFile,1}] = fileparts(deltaDataFiles{iFile});
    i             = regexp(deltaFiles{iFile}, '\d\d\d\d\d\d');
    subjID(iFile) = str2double(deltaFiles{iFile}(i:i+5));
    i2            = regexp(deltaFiles{iFile}, 'MEG_\d');
    fileID(iFile) = str2double(deltaFiles{iFile}(i+4));
    
    fprintf('Computing PLI for file %d of %d. \n', iFile, nFiles);
    
    % load in data
    tmp = load(deltaDataFiles{iFile});
    nodeData = ROInets.demean(tmp.nodeData,2);
    clear tmp;
    
    % run PLI
    [deltaPLI(:,:,iFile), deltawPLI(:,:,iFile), deltaPLV(:,:,iFile), deltaCoh(:,:,iFile), deltaIMC(:,:,iFile), deltaPSI(:,:,iFile), deltaMI(:,:,iFile), deltaPLT(:,:,iFile), deltaIMPC(:,:,iFile), deltaPCoh(:,:,iFile), deltaDTF(:,:,iFile), deltaPDC(:,:,iFile)] = PLI_band2(nodeData, Fs, windowLength, [1 4], leakageCorrected);
end%for
save(fullfile(resultsDir, 'deltaPLImetrics'), 'deltaPLI', 'deltawPLI', 'deltaPLV', 'deltaCoh', 'deltaIMC', 'deltaPSI', 'deltaMI', 'deltaPLT', 'deltaIMPC', 'deltaPCoh', 'deltaDTF', 'deltaPDC', 'subjID', 'fileID');

for iFile = length(thetaDataFiles):-1:1,
    % get file info
    [~,thetaFiles{iFile,1}] = fileparts(thetaDataFiles{iFile});
    i             = regexp(thetaFiles{iFile}, '\d\d\d\d\d\d');
    subjID(iFile) = str2double(thetaFiles{iFile}(i:i+5));
    i2            = regexp(thetaFiles{iFile}, 'MEG\d');
    fileID(iFile) = str2double(thetaFiles{iFile}(i+4));
    
    fprintf('Computing PLI for file %d of %d. \n', iFile, nFiles);
    % load in data
    tmp = load(thetaDataFiles{iFile});
    nodeData = tmp.nodeData;
    clear tmp;
    % run PLI
    [thetaPLI(:,:,iFile), thetawPLI(:,:,iFile), thetaPLV(:,:,iFile), thetaCoh(:,:,iFile), thetaIMC(:,:,iFile), thetaPSI(:,:,iFile), thetaMI(:,:,iFile), thetaPLT(:,:,iFile), thetaIMPC(:,:,iFile), thetaPCoh(:,:,iFile), thetaDTF(:,:,iFile), thetaPDC(:,:,iFile)] = PLI_band2(nodeData, Fs, windowLength, [4 8], leakageCorrected);
end%for
save(fullfile(resultsDir, 'thetaPLImetrics'), 'thetaPLI', 'thetawPLI', 'thetaPLV', 'thetaCoh', 'thetaIMC', 'thetaPSI', 'thetaMI', 'thetaPLT', 'thetaIMPC', 'thetaPCoh', 'thetaDTF', 'thetaPDC', 'subjID', 'fileID');
%{
for iFile = length(alphaDataFiles):-1:1,
    % get file info
    [~,alphaFiles{iFile,1}] = fileparts(alphaDataFiles{iFile});
    i             = regexp(alphaFiles{iFile}, '\d\d\d\d\d\d');
    subjID(iFile) = str2double(alphaFiles{iFile}(i:i+5));
    i2            = regexp(alphaFiles{iFile}, 'MEG\d');
    fileID(iFile) = str2double(alphaFiles{iFile}(i+4));
    
    fprintf('Computing PLI for file %d of %d. \n', iFile, nFiles);
    % load in data
    tmp = load(alphaDataFiles{iFile});
    nodeData = tmp.nodeData;
    clear tmp;
    % run PLI
    [alphaPLI(:,:,iFile), alphawPLI(:,:,iFile), alphaPLV(:,:,iFile), alphaCoh(:,:,iFile), alphaIMC(:,:,iFile), alphaPSI(:,:,iFile), alphaMI(:,:,iFile), alphaPLT(:,:,iFile), alphaIMPC(:,:,iFile), alphaPCoh(:,:,iFile), alphaDTF(:,:,iFile), alphaPDC(:,:,iFile)] = PLI_band2(nodeData, Fs, windowLength, [8 13], leakageCorrected);
end%for
save(fullfile(resultsDir, 'alphaPLImetrics'), 'alphaPLI', 'alphawPLI', 'alphaPLV', 'alphaCoh', 'alphaIMC', 'alphaPSI', 'alphaMI', 'alphaPLT', 'alphaIMPC', 'alphaPCoh', 'alphaDTF', 'alphaPDC', 'subjID', 'fileID');
%}


for iFile = length(betaDataFiles):-1:1,
    % get file info
    [~,betaFiles{iFile,1}] = fileparts(betaDataFiles{iFile});
    i             = regexp(betaFiles{iFile}, '\d\d\d\d\d\d');
    subjID(iFile) = str2double(betaFiles{iFile}(i:i+5));
    i2            = regexp(betaFiles{iFile}, 'MEG\d');
    fileID(iFile) = str2double(betaFiles{iFile}(i+4));
    
    fprintf('Computing PLI for file %d of %d. \n', iFile, nFiles);
    % load in data
    tmp = load(betaDataFiles{iFile});
    nodeData = tmp.nodeData;
    clear tmp;
    % run PLI
    [betaPLI(:,:,iFile), betawPLI(:,:,iFile), betaPLV(:,:,iFile), betaCoh(:,:,iFile), betaIMC(:,:,iFile), betaPSI(:,:,iFile), betaMI(:,:,iFile), betaPLT(:,:,iFile), betaIMPC(:,:,iFile), betaPCoh(:,:,iFile), betaDTF(:,:,iFile), betaPDC(:,:,iFile)] = PLI_band2(nodeData, Fs, windowLength, [13 30], leakageCorrected);
end
save(fullfile(resultsDir, 'betaPLImetrics'),  'betaPLI', 'betawPLI', 'betaPLV', 'betaCoh', 'betaIMC', 'betaPSI', 'betaMI', 'betaPLT', 'betaIMPC', 'betaPCoh', 'betaDTF', 'betaPDC', 'subjID', 'fileID');


for iFile = length(gammaDataFiles):-1:1,
    % get file info
    [~,gammaFiles{iFile,1}] = fileparts(gammaDataFiles{iFile});
    i             = regexp(gammaFiles{iFile}, '\d\d\d\d\d\d');
    subjID(iFile) = str2double(gammaFiles{iFile}(i:i+5));
    i2            = regexp(gammaFiles{iFile}, 'MEG\d');
    fileID(iFile) = str2double(gammaFiles{iFile}(i+4));
    
    fprintf('Computing PLI for file %d of %d. \n', iFile, nFiles);
    % load in data
    tmp = load(gammaDataFiles{iFile});
    nodeData = tmp.nodeData;
    clear tmp;
    % run PLI
    [gammaPLI(:,:,iFile), gammawPLI(:,:,iFile), gammaPLV(:,:,iFile), gammaCoh(:,:,iFile), gammaIMC(:,:,iFile), gammaPSI(:,:,iFile), gammaMI(:,:,iFile), gammaPLT(:,:,iFile), gammaIMPC(:,:,iFile), gammaPCoh(:,:,iFile), gammaDTF(:,:,iFile), gammaPDC(:,:,iFile)] = PLI_band2(nodeData, Fs, windowLength, [13 30], leakageCorrected);
end
save(fullfile(resultsDir, 'gammaPLImetrics'),  'gammaPLI', 'gammawPLI', 'gammaPLV', 'gammaCoh', 'gammaIMC', 'gammaPSI', 'gammaMI', 'gammaPLT', 'gammaIMPC', 'gammaPCoh', 'gammaDTF', 'gammaPDC', 'subjID', 'fileID');
% save(fullfile(resultsDir, 'PLIfiles'), 'gammaDataFiles', '-append');
save(fullfile(resultsDir, 'PLIfiles'), 'alphaDataFiles', 'betaDataFiles', 'thetaDataFiles', 'deltaDataFiles', 'gammaDataFiles');


