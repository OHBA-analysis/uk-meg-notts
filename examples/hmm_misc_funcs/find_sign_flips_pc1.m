function [ Dsnew, flips ] = find_sign_flips_pc1( Ds )

% netmats
num_subj=length(Ds);

Dp = spm_eeg_load(Ds{1});
num_nodes=size(Dp,1);

flips=ones(num_nodes,num_subj);

for subnum = 1:num_subj       
    disp(['Computing for subj num ' num2str(subnum)]);
       
    Dp = spm_eeg_load(Ds{subnum});
                    
    % no embedding
    embed.do=0;
    embed.centre_freq=13;
    embed.tres=1/Dp.fsample;
    roinets_protocol='none';
    normalisation='voxelwise';
    logtrans=0;

    % returns data as num_nodes x num_embeddings x ntpts
    datap = prepare_data(Dp,normalisation,logtrans,embed,roinets_protocol);        
    datap=permute(datap,[1 3 2]);
      
    % compute pc1
    pcadim = 1;
    [allsvd,M] = eigdec(cov(datap'),pcadim);
    pc1 = M' * datap; % should be 1 x ntpts
    
    for nn=1:num_nodes,
        ccs=corrcoef(datap(nn,:),pc1);
        if ccs(1,2)<0
            flips(nn,subnum)=-1;
        end;
    end;

    % create new object
    Dsnew{subnum}=Dp;
    S=[];
    S.D=Dsnew{subnum};
    S.outfile = prefix(Ds{subnum},'flip');
    Dsnew{subnum}=spm_eeg_copy(S);
    
    for pp=1:size(Dp,1)
        for tri=1:size(Dsnew{subnum},3)
            Dsnew{subnum}(pp,:,tri)=Dp(pp,:,tri)*flips(pp,subnum);            
        end;       
    end;
    
    Dsnew{subnum}.save;
end;
