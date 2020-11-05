function [ res ] = synchro_hmm_multichan( S )

% [ res ] = synchro_hmm_multichan( S )
%
% MWW 2014

if ~isfield(S,'norm_vectors')
    S.norm_vectors=0;
end;

if ~isfield(S,'orth_node_tcs')
    S.orth_node_tcs=0;
end;

tres=S.tres;

fsample=1/tres;

%deltasecs=0.0025;
deltasecs=S.deltasecs;
delta=(deltasecs/tres);
delta=ceil(delta);

Nchans=size(S.data,2);

if size(S.data,2)>size(S.data,1),
    error('Data should be num_nodes x num_tpts');
end;

xin=normalise(S.data)';

if(S.orth_node_tcs && size(xin,1)>1)

    xin = get_orthogonalised_node_tcs(xin, eye(Nchans), 'simpleOrthPCs','symmetric');
end;

xin=normalise(xin')';

M=S.M;

dim1=size(xin(:,delta*(M-1)+1:end),1);
dim2=size(xin(:,delta*(M-1)+1:end),2);
xe=zeros(M*dim1,dim2);

dim1=1;
to=0;
for chan=1:Nchans,
    fromstart=to+1;
    for m=0:M-1,
        %m/(M-1)
        from=m*dim1+fromstart;
        to=(fromstart-1)+m*dim1+dim1;

        xe(from:to,:)=xin(chan,delta*(M-1-m)+1:end-delta*m);    
    end;
end;

if S.norm_vectors>0,
    for ii=1:size(xe,2),
        
        if S.norm_vectors>1,
            [tmp ai]=sort(xe(:,ii));
            xe(:,ii)=ai;
        end;
        
        xe(:,ii)=normalise(xe(:,ii));
        
    end;
end;

S2=[];
S2.data=[xe]';

if(S.hmm_pca_dim<=0)
    S.hmm_pca_dim=size(S2.data,2);
end;

%% do PCA
disp('doing pre-HMM PCA');
[allsvd,Apca]=pca(S2.data,S.hmm_pca_dim);
pinvApca=pinv(Apca);
res.Apca=Apca;

S2.data=(pinvApca*S2.data')';   

S2.data=normalise(S2.data);
S2.NK=S.NK;
S2.force_zero_means=S.force_zero_means;
S2.num_starts=S.num_hmm_starts;
disp('doing HMM');
[ hmm, block, frenbest ] = run_multistart_hmm( S2 );

res.xin=xin;
res.hmm=hmm;
res.block=block;
res.delta=delta;
res.S=S;
res.frenbest=frenbest;

if(S.do_plots),
    
    NK=S.NK;
    figure; hold on;
    for ii=1:NK,
        %ts=tres:tres:length(block(1).q_star)*tres;
        ts=tres:tres:size(res.xin,2)*tres;        
        ts=ts(delta*(M-1)+1:end);
        plot(ts,(res.block(1).q_star==ii)+2*(ii) ,'b');
        xlabel('time(s)');
    end;
    title(S.title);

    figure;
    for ii=1:NK,

        cv=res.Apca*res.hmm.state(ii).Cov*res.Apca';
        cv=corrcov(cv);
        subplot(NK,1,ii); imagesc(abs(cv),[0 1]); colorbar; 
%        subplot(NK,1,(ii-1)*2+2); imagesc(abs(res.hmm.state(ii).Mu), [0 1]); colorbar; 

    end;
    title(S.title);

    
    StatePath = res.block(1).q_star;
    %StatePath = simdata.Xclass(1:Nto);
    NumberOccurences = zeros(1,res.hmm.K);
    FractionalOccupancy = zeros(1,res.hmm.K);
    MeanLifeTime = zeros(1,res.hmm.K);

    for ctr = 1:res.hmm.K
        temp = StatePath == ctr;
        NumberOccurences(ctr) = sum(diff(temp) == 1);
        FractionalOccupancy(ctr) = sum(temp)/length(temp); 
        MeanLifeTime(ctr) = sum(temp)/NumberOccurences(ctr);
    end
    x=1:NK;

    % plots
    figure;
    subplot(121);bar(x,FractionalOccupancy);
    plot4paper('State #','Fractional occupancy');

    subplot(122);bar(x,(tres)*MeanLifeTime);
    plot4paper('State #','Mean life time (s)');
    title(S.title);

end;

end

