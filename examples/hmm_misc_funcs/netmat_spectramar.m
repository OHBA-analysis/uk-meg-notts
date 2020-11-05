function [ state ] = netmat_spectramar( data, S )

% data is num_nodes x num_embeddings x ntpts
% netmat is num_nodes x num_nodes x num_embeddings
% netmat_full is (num_nodes x num_embeddings) x (num_nodes x num_embeddings)
   
try S.type=S.type; catch S.type='coh'; end
try S.var_normalise=S.var_normalise; catch S.var_normalise=false; end;

global OSLDIR;


if isfield(S,'netmats') && isfield(S.netmats,'spectramt') && ~isempty(S.netmats.spectramt)
    netmats=S.netmats;
else

    if size(data,2)>1
        error('spectramt only compatible with 1 embedding (i.e. the raw time series)');
    end;
    
    Hz=S.fsample;
    
    if S.var_normalise
        data=normalise(data,3);
    end;

    data=permute(data,[3 1 2]);
    
    
    % compute spectra
    params = struct('Fs',Hz); % Sampling rate
    params.fpass = S.fband;  % band of frequency you're interested in 
         
    params.completelags = 1; % complete lags? default to 1 (recommended)
    params.MLestimation = 1; % maximum likelihood or Bayesian - former is recommended
    params.order=S.order;
    fitmt = hmmspectramar(data,Ts,hmm,gamma,params);

    for kk=1:length(fitmt.state)
        state{kk}.spectramt=fitmt.state(kk);
    end
    
end

for kk=1:length(state)

    % by default set netmat to coh, with num_embeddings in output
    % corresponding to num freq bins: d
    tmp=state{kk}.spectramt.(S.type);
    state{kk}.netmat=permute(tmp,[2 3 1]);

    % create netmat_full:

    switch S.full_type
        case 'full' % netmat_full is (num_nodes x num_embeddings) x (num_nodes x num_embeddings)
            num_embeddings=size(state{kk}.netmat,3);
            num_nodes=size(state{kk}.netmat,1);
            state{kk}.netmat_full=speye(num_nodes*num_embeddings);
            for node_ind=1:num_nodes
                from = (node_ind-1)*num_embeddings+1;
                to = from+num_embeddings-1;
                for node_ind2=1:num_nodes
                    from2 = (node_ind2-1)*num_embeddings+1;
                    to2 = from2+num_embeddings-1;
                    state{kk}.netmat_full(from:to,from2:to2)=sparse(diag(permute(state{kk}.netmat(node_ind,node_ind2,:),[3 1 2])));
                end;
            end;
        case 'mean_abs' % netmat_full is (num_nodes x num_nodes)

            state{kk}.netmat_full=mean(abs(state{kk}.netmat),3);

    end
end

netmats.type=S.type;
netmats.full_type=S.full_type;

end

