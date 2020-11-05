%% main run hmm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP THE MATLAB PATHS
% make sure that fieldtrip and spm are not in your matlab path
    
tilde='/home/mwoolrich/';
tilde='/Users/woolrich/';

cd([tilde '/homedir/scripts/osl/osl-core']);
%cd /Users/woolrich/homedir/osl/current_release/osl/osl-core

osl_startup();

addpath(genpath([tilde 'Dropbox/vols_scripts/hmm_misc_funcs']));
addpath(genpath([tilde 'Dropbox/vols_scripts/notts_ukmp/main']));

% sessionname='eo'; nohup /Applications/MATLAB_R2018b.app/bin/matlab -nodesktop -nodisplay < /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/main/main_run_hmm.m > /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/logs/main_run_hmm_$sessionname.log 2>&1 &
% sessionname='eo'; tail -f /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/logs/main_run_hmm_$sessionname.log

% more /Users/woolrich/Dropbox/vols_scripts/notts_ukmp/main/main_run_hmm.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

spmfilesdir='/Users/woolrich/homedir/vols_data/notts_ukmp/spm/';

preproc_name='sept2019'; 
session_name='eo';nsubjects=55;

%session_name='vml';nsubjects=51;
%session_name='vml';nsubjects=2;
%session_name='vms';nsubjects=52;


prep_sessions_to_do=1:nsubjects;
sessions_to_exclude=[];
hmm_sessions_to_do=setdiff(prep_sessions_to_do,sessions_to_exclude);

%%

S=[];

S.do_prepare=1; 
S.do_hmm=1; 
S.do_spectral_estimation=1;
S.preproc_name=preproc_name;
S.session_name=session_name; 
S.spmfilesdir=spmfilesdir;
S.prep_sessions_to_do=prep_sessions_to_do;
S.hmm_sessions_to_do=hmm_sessions_to_do;
S.num_embeddings=13;
S.freq_range=[1 45];

S.hmm_name='';
disp(S);

[hmm hmmfname hmmoptions settings_prepare] = run_full_hmm_pipeline(S);

disp(hmmfname);
