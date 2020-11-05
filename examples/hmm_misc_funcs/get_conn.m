function conn=get_conn(netmat_full,netmat,node_ind,node_ind2,num_embeddings)

% conn=get_conn(netmat_full,netmat,node_ind,node_ind2,num_embeddings)
% returns num_embeddings x num_embeddings

from = (node_ind-1)*num_embeddings+1;
to = from+num_embeddings-1;

from2 = (node_ind2-1)*num_embeddings+1;
to2 = from2+num_embeddings-1;

conn=netmat_full(from:to,from2:to2);

end

