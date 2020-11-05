%% main_preproc
%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path
    
tilde='/home/mwoolrich/';
tilde='/Users/woolrich/';

cd([tilde '/homedir/scripts/osl/osl-core']);
%cd /Users/woolrich/homedir/osl/current_release/osl/osl-core

osl_startup();

addpath(genpath([tilde 'Dropbox/vols_scripts/hmm_misc_funcs']));
addpath(genpath([tilde 'Dropbox/vols_scripts/notts_ukmp/main']));

% MRI conversions need freesurfer:
addpath('/Applications/freesurfer/matlab/');

% nohup /Applications/MATLAB_R2017b.app/bin/matlab -nodesktop -nodisplay < /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/main/main_preproc.m > /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/logs/main_preproc_eo.log 2>&1 &
% tail -f /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/logs/main_preproc_eo.log

% nohup /Applications/MATLAB_R2017b.app/bin/matlab -nodesktop -nodisplay < /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/main/main_preproc.m > /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/logs/main_preproc_vml.log 2>&1 &
% tail -f /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/logs/main_preproc_vml.log

%%%%%%%%%%%%%%%%%%
%% Settings

do_convert=true;

do_mri_convert=true;

use_prerun_opt=false;

session_name='eo';
%session_name='vml';
%session_name='vms';

preproc_name='oct2019';

%%%%%%%%%%%%%%%%%%
%% SET UP FILE NAMES

workingdir  = [tilde 'homedir/vols_data/notts_ukmp/'];
spmfilesdir = fullfile(workingdir,'spm/');

datadir     = [tilde 'homedir/vols_data/notts_ukmp/raw_data/'];
africadir   = [spmfilesdir 'africa/'];

subjectdirs = dir([datadir '/3*']);
subjectdirs = sort_nat({subjectdirs.name});

ctf_files_eo            = cell(size(subjectdirs));
ctf_files_vmg           = cell(size(subjectdirs));
ctf_files_vml           = cell(size(subjectdirs));
ctf_files_vms           = cell(size(subjectdirs));
spm_files_eo            = cell(size(subjectdirs));
spm_files_vmg           = cell(size(subjectdirs));
spm_files_vml           = cell(size(subjectdirs));
spm_files_vms           = cell(size(subjectdirs));
structural_files        = cell(size(subjectdirs));
pos_files               = cell(size(subjectdirs));

for s = 1:length(subjectdirs)
   
    dsfile = dir([datadir subjectdirs{s} '/*Eyes_Open*.ds']); 
    if ~isempty(dsfile)
        ctf_files_eo{s} = [datadir subjectdirs{s} '/' dsfile.name];

        % set up a list of SPM MEEG object file names (we only have one here)
        spm_files_eo{s}    = [spmfilesdir subjectdirs{s} '_eo.mat'];
    end
    
    dsfile = dir([datadir subjectdirs{s} '/*VisMotor_Gamma*.ds']); 
    if ~isempty(dsfile)
        ctf_files_vmg{s} = [datadir subjectdirs{s} '/' dsfile.name];

        % set up a list of SPM MEEG object file names (we only have one here)
        spm_files_vmg{s}    = [spmfilesdir subjectdirs{s} '_vmg.mat'];
    end
    
    dsfile = dir([datadir subjectdirs{s} '/*VisMotor_Long*.ds']); 
    if ~isempty(dsfile)
        ctf_files_vml{s} = [datadir subjectdirs{s} '/' dsfile.name];

        % set up a list of SPM MEEG object file names (we only have one here)
        spm_files_vml{s}    = [spmfilesdir subjectdirs{s} '_vml.mat'];
    end
 
    
    dsfile = dir([datadir subjectdirs{s} '/*VisMotor_Short*.ds']); 
    if ~isempty(dsfile)
        ctf_files_vms{s} = [datadir subjectdirs{s} '/' dsfile.name];

        % set up a list of SPM MEEG object file names (we only have one here)
        spm_files_vms{s}    = [spmfilesdir subjectdirs{s} '_vms.mat'];
    end
        
    % structural files
    try
        
        niifile = [datadir subjectdirs{s} '/' subjectdirs{s} '_CRG']; 
        
        if do_mri_convert || isempty(dir([niifile '.nii']))
            runcmd(['rm -f ' niifile '.nii']);
            runcmd(['rm -f ' niifile '.nii.gz']);
            runcmd(['rm -f ' datadir subjectdirs{s} '/y_' subjectdirs{s} '_CRG.nii']);
            runcmd(['rm -f ' niifile '_*.gii']);
            runcmd(['rm -f ' niifile '*_*.*']);
            
            mrifile = dir([datadir subjectdirs{s} '/' num2str(subjectdirs{s}) '_CRG.mri']); 
            
            if isempty(mrifile)
                me=MException('Preprocess:noStuctural',['No mri or nii structural file for subject ' subjectdirs{s}]);
                throw(me);
            end           
            
            %[ niifile ] = mri2analyze( [datadir subjectdirs{s} '/' mrifile.name] );
            
            
            niifilename=[niifile '.nii'];
            convert_mri([datadir subjectdirs{s} '/' mrifile.name], niifilename);
           
        end
        
        niifilename=[niifile '.nii'];
        structural_files{s} = niifilename;
        
    catch me
        warning(me.message)
        structural_files{s} = [];
    end
    
    % list of head position files
    
    pfname=[datadir 'pos_files/' subjectdirs{s} '.pos'];
    pf = dir(pfname); 
    if ~isempty(pf)
        pos_files{s}=pfname;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Select spm files to work with

switch session_name
    case 'eo'
        spm_files=spm_files_eo;
        ctf_files=ctf_files_eo;
        event_type=[];
    case 'vml'
        spm_files=spm_files_vml;
        ctf_files=ctf_files_vml;
        event_type={'Abduction','Stim_Onset'};
    case 'vms'
        spm_files=spm_files_vms;
        ctf_files=ctf_files_vms;
        event_type={'Abduction','Stim_Onset'};
    case 'vmg'
        spm_files=spm_files_vmg;
        ctf_files=ctf_files_vmg;
        event_type={'Abduction','Stim_Onset'};
end

subjects_to_do = find(~cellfun(@isempty,spm_files) & ~cellfun(@isempty,structural_files));
subjects_to_do=subjects_to_do(1);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% CTF data preprocessing
%% CONVERT FROM .ds TO AN SPM MEEG OBJECT:

if do_convert

    S2=[];
    mkdir([tilde '/homedir/vols_data/notts_ukmp/results/']);

    for ss=1:length(subjects_to_do), % iterates over subjects    
        s=subjects_to_do(ss);

        S2.outfile=spm_files{s};

        % The conversion to SPM will show a histogram of the event codes
        % and correspond to those listed below in the epoching section
        D = osl_import(ctf_files{s},S2);

    end

    %%%%%%%%%%%%%%%%%%%
    % Sort out fiducials:

    for ss=1:length(subjects_to_do) % iterates over subjects    
        f=subjects_to_do(ss);

        spm_file = prefix(spm_files{f},'');
        D = spm_eeg_load(spm_file);

        fID = fopen(pos_files{f});
        fid_data = textscan(fID, '%s %f %f %f');
        fclose(fID);

        fid_new = [];

        % Fiducials:
        fid_inds = [find(strcmpi(fid_data{1},'nasion'))
            find(strcmpi(fid_data{1},'left'))
            find(strcmpi(fid_data{1},'right'))];

        fid_new.fid.pnt = [fid_data{2}(fid_inds) fid_data{3}(fid_inds) fid_data{4}(fid_inds)] * 10;
        fid_new.fid.label = {'nas';'lpa';'rpa'};

        % Headshape:
        hs_inds = setdiff(2:length(fid_data{1}),fid_inds);
        fid_new.pnt = [fid_data{2}(hs_inds) fid_data{3}(hs_inds) fid_data{4}(hs_inds)] * 10;
        %nose = fid_new.pnt(:,3) < -10;
        %fid_new.pnt(nose,:) = [];

        fid_new.unit = 'mm';

        % Labels:
        fid_labels = fid_data{1};

        D = fiducials(D,fid_new);
        D.save;

    end
end

%%%%%%%%%%%%%%%%%%%
%% Check headshape points:
if false
    for ss=1:length(subjects_to_do) % iterates over subjects    
        f=subjects_to_do(ss);
        D = spm_eeg_load(spm_files{f});
        osl_edit_fid(D)
        waitfor(gcf)
    end
end

%%%%%%%%%%%%%%%%%%%
%% OPT

opt=[];
opt.dirname=[workingdir preproc_name '_' session_name '.opt/'];

if use_prerun_opt
    
    opt=osl_load_opt(opt.dirname);
    %opt=rmfield(opt,'osl2_version');
    opt=opt_consolidate_results(opt);

else
    
    runcmd(['rm -rf ' opt.dirname]);
    
    clear spm_files_in structural_files_in;

    for ss=1:length(subjects_to_do), % iterates over subjects  
        
        f=subjects_to_do(ss);
        
        spm_files_in{ss} = prefix(spm_files{f},'');
        structural_files_in{ss} = structural_files{f};     
        
        % label artefact chans:
        D=spm_eeg_load(spm_files_in{ss});
        D = D.chantype(find(strcmp(D.chanlabels,'EEG060')),'EMG');
        D = D.chantype(find(strcmp(D.chanlabels,'EEG059')),'ECG');
        D = D.chantype(find(strcmp(D.chanlabels,'EEG057')),'EOG1');
        D = D.chantype(find(strcmp(D.chanlabels,'EEG058')),'EOG2');
        D.save;

    end

    opt.datatype='ctf';
    opt.spm_files=spm_files_in;

    opt.downsample.do=1;
    opt.downsample.freq=250;
    
    opt.highpass.do = 1;
    opt.highpass.cutoff = 0.1;
    
    opt.mains.do = 0;
          
    opt.africa.todo.ica=0;
    opt.africa.todo.ident=0;
    opt.africa.todo.remove=0;
    
    opt.africa.ident.func = @identify_artefactual_components_auto;
    opt.africa.ident.kurtosis_wthresh=0.2;
    opt.africa.ident.max_num_artefact_comps=2;
    opt.africa.precompute_topos=1;

    opt.bad_segments.do=~(length(event_type)>0);
    opt.bad_segments.event_significance=0.05;

    opt.coreg.do=true;
    opt.coreg.mri=structural_files_in;
    opt.coreg.use_rhino=1;
    opt.coreg.useheadshape=1;
    opt.coreg.forward_meg='MEG Local Spheres';
    
    % Epoching settings
    %
    % Here the epochs are set to be from -1s to +2s relative to triggers
    % in the MEG data.
    opt.epoch.do=(length(event_type)>0);
    if opt.epoch.do
        for ee=1:length(event_type)
            opt.epoch.trialdef(ee).conditionlabel = event_type{ee};
            opt.epoch.trialdef(ee).eventtype = event_type{ee};
            opt.epoch.trialdef(ee).eventvalue = 1;
        end
        
        switch session_name
            case 'vms'
                opt.epoch.time_range=[-4 4];
            otherwise
                opt.epoch.time_range=[-10 10];
        end
        opt.outliers.do = true;
        opt.outliers.outlier_measure_fns = {'std'};
        %opt.outliers.event_significance=0.05;
    else
        opt.outliers.do = 0;
    end
   
    % opt=opt_consolidate_results(opt);
    opt=osl_run_opt(opt);
end

if 0

    %%%%%%%%%%%%%%%%%%%
    %% Check & correct coregistration:
    for ss=1:length(subjects_to_do), % iterates over subjects    
        f=subjects_to_do(ss);
        %D = spm_eeg_load(prefix(spm_files{f},'A'));
        D = spm_eeg_load(opt.results.spm_files{ss});
        if S.use_rhino
            rhino_manual(D);
        else
            spm_eeg_inv_checkdatareg(D);
            %spm_eeg_inv_checkmeshes(D, 1);
        end
        waitfor(gcf)
    end

    %%%%%%%%%%%%%%%%%%%
    %% Rerun forward model if position edited in rhino_manual
    for ss=1:length(subjects_to_do), % iterates over subjects    
        f=subjects_to_do(ss);
        S.D             = opt.results.spm_files{ss}
        S.forward_meg   = 'MEG Local Spheres';
        osl_forward_model(S);
    end
end

%%%%%%%
%% copy and rename spm files to have names session1...

% session 8 gets its own state!

switch session_name
    case 'eo'
        sessions_to_do=[1:7, 9:length(opt.results.spm_files)];
    case 'vml'
        sessions_to_do=[1:7, 9:length(opt.results.spm_files)];
        sessions_to_do=1:2;
    case 'vms'
        sessions_to_do=[1:7, 9:length(opt.results.spm_files)];
        sessions_to_do=1:3;
end

spm_links=cell(length(sessions_to_do),2);
for ss=1:length(sessions_to_do),

    S=[];
    S.D=opt.results.spm_files{sessions_to_do(ss)};
    S.outfile=[spmfilesdir preproc_name '_' session_name '_session' num2str(ss)];
    Dnew=spm_eeg_copy(S);
    
    % write log of how original spm files are linked to copies
    spm_links{ss,1}=S.D;
    spm_links{ss,2}=S.outfile;
    
    % if epoched data:
    switch session_name
        case {'vml','vms'}
    
            S=[];
            S.D=opt.results.spm_files_epoched{sessions_to_do(ss)};
            S.outfile=[spmfilesdir preproc_name '_' session_name '_session' num2str(ss) '_epoched'];
            Dnew=spm_eeg_copy(S);

    end

end

%save([spmfilesdir 'spm_links_' preproc_name '_' session_name ], 'spm_links');

%%

if 0
    %%
    structural='/Users/woolrich/homedir/vols_data/notts_ukmp/raw_data/3004/3004_CRG.nii';
    Dfile='/Users/woolrich/homedir/vols_data/notts_ukmp/spm/3004_eo'
    % fsleyes /Users/woolrich/homedir/vols_data/notts_ukmp/raw_data/3004/3004_CRG.nii

    S = [];
    S.D                 = Dfile
    S.mri               = structural
    S.useheadshape      = 1;
    S.use_rhino         = 1;
    S.forward_meg       = 'MEG Local Spheres';
    S.fid.label.nasion  = 'nas';
    S.fid.label.lpa     = 'lpa';
    S.fid.label.rpa     = 'rpa';
    D = osl_headmodel(S);
end

