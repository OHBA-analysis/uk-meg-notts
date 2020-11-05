function [ state_netmats ] = hmm_state_netmats( hmm, S )

% [ state_netmats ] = hmm_state_netmats( hmm, S )

try S.global_only=S.global_only; catch S.global_only=false; end
try parcellated_filenames=S.parcellated_filenames; catch error('Need to specify parcellated_filenames'); end

logtrans=0;

normalisation  = S.concat.normalisation;
roinets_protocol=S.concat.protocol;
embed=S.concat.embed;
assignment=S.assignment;

if isempty(hmm)
    
    D_fnames=S.D_fnames;
    %parcelAssignments=S.parcelAssignments;
    
else
    
    D_fnames=S.parcellated_filenames;
    
end;

nsubs=length(D_fnames);

for subnum = 1:nsubs

    disp(['Computing for subj num ' num2str(subnum)]);
       
    try
        Dp = spm_eeg_load(prefix(D_fnames{subnum},'p'));
    catch
        Dp = spm_eeg_load((D_fnames{subnum}));
    end;
        
    %state_netmats{subnum}.parcelAssignments=parcelAssignments;
    %state_netmats{subnum}.parcellation=Dp.parcellation;
    
    embed.tres=1/Dp.fsample;
        
    %     if subnum==3,
    %         Dp(2,:,:)=-Dp(2,:,:);
    %         Dp(37,:,:)=-Dp(37,:,:);
    %         Dp.save;
    %     end;

    % returns data as num_nodes x num_embeddings x ntpts
    datap = prepare_data(Dp,normalisation,logtrans,embed,roinets_protocol);        
    
    if ~S.global_only && hmm.K>1
        
        hmm_sub = hmm; 
        hmm_sub.statepath = hmm.statepath(hmm.subj_inds==subnum); 
%        hmm_sub.statepath_soft = hmm.train.Gamma(subj_inds==subnum,:);

%        hmm_sub = rmfield(hmm_sub,'MixingMatrix');    

        for k = 1:hmm_sub.K

             disp(['Computing for state ' num2str(k)]);

             switch assignment
                 case 'hard' 
                    if 1,
                        inds = logical(hmm_sub.statepath == k);
                        if size(datap,3) ~= length(inds)
                            % need to resample state time course
                            inds=logical(round(resample(double(inds),size(datap,3),length(inds))));               
                        end;
                        normconstant = sqrt(size(inds,1) / sum(double(inds)) );
                        
                        datap_in=datap(:,:,inds) * normconstant;
                    else
                        x = double(hmm_sub.statepath == k);
                        if size(datap,3) ~= length(x)
                            % need to resample state time course
                            error('not implemented');
                        end;

                        normconstant = sqrt(size(x,1) / sum(x) ); 
                        x = x * normconstant; 
                        x2=permute(repmat(x,[1, size(datap,1), size(datap,2)]),[2 3 1]);
                        datap_in=datap.*x2;
                    end;
                     
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

            [state_netmats{subnum}.state{k}]=feval(S.netmat_method,datap_in,S.netmat_method_options);

        end;
    end;
    
    [state_netmats{subnum}.global]=feval(S.netmat_method,datap,S.netmat_method_options);
 
    %sum(state_netmats{subnum}.global.spectramt.psd(:,9,9),1)
    %sum(state_netmats{subnum}.state{k}.spectramt.psd(:,9,9),1)
    %std(squeeze(datap_in(9,1,:))),std(squeeze(datap(9,1,:)))
    %sum(state_netmats{subnum}.global.netmat(9,9,:))
    %sum(state_netmats{subnum}.state{k}.netmat(9,9,:))
    
    
    state_netmats{subnum}.netmat_method=S.netmat_method;
end                        
