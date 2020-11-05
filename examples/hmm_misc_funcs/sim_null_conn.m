function [null_conn_mean null_conn_std]=sim_null_conn(data_files, hmm, embed_centre_freq)

%% create null data files for use later
clear null_files;
for dd=1:length(data_files)
    D=spm_eeg_load(data_files{dd});
    null_files{dd}=prefix(D.fullfile,'null_');
    Dnull=D.copy(null_files{dd});
end

%% precompute stds ar etc.
clear stds4null ars
NP=1;
ars=zeros(length(null_files),size(D,1),NP);
for dd=1:length(null_files)
    dd
    D=spm_eeg_load(data_files{dd});
    for ii=1:size(D,1),
        stds4null(dd,ii)=std(squash(D(ii,:,:)));
        ts = (squash(D(ii,:,:)));
        
        %res=pacf(demean(ts(1:10000)),NP,NP);
        
        model = ar (demean(ts(1:10000)),NP,'burg');
        ars(dd,ii,:) = getpvec(model);
    end
    
end

save([prefix(null_files{1},'ars')],'ars');

%%

nnulls=1;

for nn=1:nnulls

    nn
    
    % create null data
    disp('create null data');
    for dd=1:length(null_files)
        
        dd
        Dnull=spm_eeg_load(null_files{dd});
        for ii=1:size(Dnull,1),
            %Dnull(ii,:,:)=randn([1,size(Dnull,2),size(Dnull,3)])*stds4null(dd,ii);

            y=spm_mar_gen(0,squeeze(ars(dd,ii,:))',1,size(Dnull,2)*size(Dnull,3));            
            Dnull(ii,:,:)=reshape(y,[1,size(Dnull,2),size(Dnull,3)]);
        end
        Dnull.save;
    end
    
    disp('compute conns');
    S=[];
    S.parcellated_filenames=null_files;
    %S.concat = hmmoptions.concat;
    S.protocol='symmetric';
    S.assignment='hard';
    S.global_only=false;

    S.netmat_method=@netmat_cov;
    S.embed.do=1;
    S.embed.rectify=false;
    S.embed.centre_freq=embed_centre_freq;
    S.normalisation='voxelwise';
    %S.normalisation='none';
    S.netmat_method_options.type='full';
    S.netmat_method_options.var_normalise=true;

    [ state_netmats_cov_null ] = hmm_state_netmats_teh( hmm, S );

    do_abs=1;
    state_netmats=state_netmats_cov_null;

    NK=length(state_netmats{1}.state);
    num_nodes=size(state_netmats{1}.state{1}.netmat,1);
    num_embeddings=size(state_netmats{1}.state{1}.netmat,3);
    conns_null=zeros(length(state_netmats),NK,num_nodes,num_nodes);
    for ss=1:length(state_netmats)
        for k=1:NK,    
            for ii=1:num_nodes,
                for jj=1:ii-1,
                    conn=get_conn(state_netmats{ss}.state{k}.netmat_full,state_netmats{ss}.state{k}.netmat,ii,jj,num_embeddings);        

                    % collapse over embeddings
                    if do_abs
                        tmp=abs(squash(conn));
                    else
                        tmp=squash(conn);
                    end

                    conns_null(ss,k,ii,jj)=mean(tmp);
                    % symmetric:
                    conns_null(ss,k,jj,ii)=conns_null(ss,k,ii,jj);
                end
            end    
        end
    end

    % normalise
    conns2=conns_null;
    for ss=1:length(state_netmats)
        for k=1:NK,    
            for ii=1:num_nodes,
                for jj=1:ii-1,
                    conns2(ss,k,jj,ii)=conns_null(ss,k,jj,ii)/(sqrt(conns_null(ss,k,jj,jj))*sqrt(conns_null(ss,k,ii,ii)));
                    % symmetric:
                    conns2(ss,k,jj,ii)=conns_null(ss,k,ii,jj);
                end
            end   
            
            % combine over edges to get distribution
            tmp=squash(triu(squeeze(conns2(ss,k,:,:))));
            inds=find(tmp~=0);
            tmp=tmp(inds);
            tmp=log((tmp)./(1-tmp));

            null_conn_mean(ss,k)=mean(tmp);
            null_conn_std(ss,k)=std(tmp);
        end
    end

    conns_null=conns2;
    
end
