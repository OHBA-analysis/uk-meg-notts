function [ state_netmats ] = hmm_state_netmats_concat( hmm, S )

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

    
end

datap=permute(datap,[1 3 2]);

if ~S.global_only && hmm.K>1
        
    hmm_sub = hmm; 
    hmm_sub.statepath = hmm.statepath(subj_inds<=nsubs); 
%        hmm_sub.statepath_soft = hmm.train.Gamma(subj_inds==subnum,:);

    hmm_sub = rmfield(hmm_sub,'MixingMatrix');    

    for k = 1:hmm_sub.K

        disp(['Computing for state ' num2str(k)]);

        switch assignment
        case 'hard' 

            inds = logical(hmm_sub.statepath == k);
            if size(datap,3) ~= length(inds)
                % need to resample state time course
                inds=logical(round(resample(double(inds),size(datap,3),length(inds))));               
            end;

            datap_in=datap(:,:,inds);

            %x = double(hmm_sub.statepath == k)';                
         case 'soft'
            x = hmm_sub.statepath_soft(:,k);

            if size(datap,3) ~= length(x)
                % need to resample state time course
                error('not implemented');
            end;

            normconstant = sqrt(size(x,1) / sum(x) ); 
            x = x * normconstant; 
            x2=permute(repmat(x,[1, size(datap,1), size(datap,2)]),[2 3 1]);
            datap_in=datap.*x2;
        end

        [state_netmats{1}.state{k}]=feval(S.netmat_method,datap_in,S.netmat_method_options);

    end;

end;
    
[state_netmats{1}.global]=feval(S.netmat_method,datap,S.netmat_method_options);

%sum(state_netmats{subnum}.global.spectramt.psd(:,9,9),1)
%sum(state_netmats{subnum}.state{k}.spectramt.psd(:,9,9),1)
%std(squeeze(datap_in(9,1,:))),std(squeeze(datap(9,1,:)))
%sum(state_netmats{subnum}.global.netmat(9,9,:))
%sum(state_netmats{subnum}.state{k}.netmat(9,9,:))


state_netmats{1}.netmat_method=S.netmat_method;
                        
