function [ netmats ] = netmat_slidingwin( data, S )

% data is num_nodes x num_embeddings x ntpts
% netmat is num_nodes x num_nodes x num_embeddings
% netmat_full is (num_nodes x num_embeddings) x (num_nodes x num_embeddings)
   
try S.type=S.type; catch S.type='coh'; end
try S.var_normalise=S.var_normalise; catch S.var_normalise=false; end;

if isfield(S,'netmats') && isfield(S.netmats,'spectramt') && ~isempty(S.netmats.spectramt)
    netmats=S.netmats;
else

    if size(data,2)>1
        error('spectramt only compatible with 1 embedding (i.e. the raw time series)');
    end;
    
    Fs=S.fsample;
    
    if S.var_normalise
        data=normalise(data,3);
    end;

    data=permute(data,[3 1 2]);

    nnodes=size(data,2);

    [netmats.spectramt.plv, netmats.spectramt.psd] = PLI_band2(data', Fs, S.windowLength, S.freqBands, false, false);

    netmats.spectramt.f=cellfun(@mean,S.freqBands);
    
end;

% by default set netmat to coh, with num_embeddings in output
% corresponding to num freq bins:
netmats.netmat=permute(netmats.spectramt.(S.type),[2 3 1]);

% create netmat_full:

switch S.full_type
    case 'full' % netmat_full is (num_nodes x num_embeddings) x (num_nodes x num_embeddings)
        num_embeddings=size(netmats.netmat,3);
        num_nodes=size(netmats.netmat,1);
        netmats.netmat_full=speye(num_nodes*num_embeddings);
        for node_ind=1:num_nodes
            from = (node_ind-1)*num_embeddings+1;
            to = from+num_embeddings-1;
            for node_ind2=1:num_nodes
                from2 = (node_ind2-1)*num_embeddings+1;
                to2 = from2+num_embeddings-1;
                netmats.netmat_full(from:to,from2:to2)=sparse(diag(permute(netmats.netmat(node_ind,node_ind2,:),[3 1 2])));
            end;
        end;
    case 'mean_abs' % netmat_full is (num_nodes x num_nodes)

        netmats.netmat_full=mean(abs(netmats.netmat),3);
        
end;

netmats.type=S.type;
netmats.full_type=S.full_type;

end

