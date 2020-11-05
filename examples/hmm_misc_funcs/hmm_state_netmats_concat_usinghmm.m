function [ state_netmats ] = hmm_state_netmats_concat_usinghmm( hmm, S )

% [ state_netmats ] = hmm_state_netmats( hmm, S )

try 
    S.global_only=S.global_only; 
catch
    S.global_only=false;
end

logtrans=0;
state_netmats=[];

normalisation  = S.concat.normalisation;
roinets_protocol=S.concat.protocol;
embed=S.concat.embed;
assignment=S.assignment;

if isempty(hmm)
    
    D_fnames=S.D_fnames;
    
else
    
    filenames=hmm.filenames;

    load(filenames.hmm);
    load(filenames.concat);
    
    D_fnames=filenames.prepare;
    
    % get parcel assignments:
    pathstr = fileparts([filenames.prepare{1}]);
    fname = fullfile(pathstr,'parcellation');
    load(fname);

end;

nsubs=length(D_fnames);

datap=[];

Ts = zeros(nsubs,1);

for subnum = 1:nsubs

    disp(['Concatenating for subj num ' num2str(subnum)]);
    try 
        D = spm_eeg_load(prefix(filenames.prepare{subnum},'p'));
    catch
        warning('Loading in SPM MEG object that might not be parcellated.');
        D = spm_eeg_load(filenames.prepare{subnum});
    end;
    
    if subnum==1,
        state_netmats{1}.parcellation=D.parcellation;
    end;

    embed.tres=1/D.fsample;

    %data = prepare_data(D,normalisation,logtrans,embed,'none');
    data = prepare_data(D,normalisation,logtrans,embed,roinets_protocol);

    data=permute(data,[1 3 2]);
    
    datap = [datap, data];

    Ts(subnum)= length(data);
end

datap=permute(datap,[1 3 2]);

if ~S.global_only && hmm.K>1
        
    % state_netmats_mt{ss}.state{k}.spectramt
    [state_netmats{1}]=feval(S.netmat_method,datap,hmm,Ts,S.netmat_method_options);

end
   
% set all to one state:
hmm_global=hmm;
hmm_global.train.Gamma=ones(size(hmm_global.train.Gamma,1),1);
hmm_global.K=1;
hmm_global.train.K=1;
hmm_global.state=hmm.state(1:1);
[tmp]=feval(S.netmat_method,datap,hmm_global,Ts,S.netmat_method_options);
state_netmats{1}.global=tmp{1};

%sum(state_netmats{subnum}.global.spectramt.psd(:,9,9),1)
%sum(state_netmats{subnum}.state{k}.spectramt.psd(:,9,9),1)
%std(squeeze(datap_in(9,1,:))),std(squeeze(datap(9,1,:)))
%sum(state_netmats{subnum}.global.netmat(9,9,:))
%sum(state_netmats{subnum}.state{k}.netmat(9,9,:))


state_netmats{1}.netmat_method=S.netmat_method;
                        
