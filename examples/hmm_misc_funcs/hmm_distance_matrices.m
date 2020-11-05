function [ distance_matrices ] = hmm_distance_matrices( S )

% [ distance_matrices ] = hmm_distance_matrices( S )
%
% Computes state x state matrix for each node indicating
% connectivity profile distances between all states
%
% state_netmats.netmat{k} is num_nodes x num_nodes x num_embeddings
% state_netmats.netmat_full{k} is (num_nodes x num_embeddings) x (num_nodes x num_embeddings)
%       where the (num_nodes x num_embeddings) dimension is in the order:
%       [node1,embedding1; node1,embedding2 ... node1,embeddingN;
%       [node2,embedding1; node2,embedding2 ... node2,embeddingN;
%       etc...
%
% nodemap{k} is num_nodes x 1 

state_netmats=S.state_netmats;

num_embeddings=size(state_netmats.netmat_global,3);
num_nodes=size(state_netmats.netmat_global,1);

global_conn_profiles=reshape(state_netmats.netmat_global,[num_nodes, (num_nodes*num_embeddings)]);
    
NK=length(state_netmats.netmat);

remove_auto=1;

clear distance_matrices;
for jj=1:num_nodes,
        
    for kk=1:NK,
        % extract state profile for this node                 
        conn_profilekk=squash(get_conn_profile(state_netmats.netmat_full{kk},jj,num_embeddings,remove_auto));

        for ii=1:NK,
            conn_profileii=squash(get_conn_profile(state_netmats.netmat_full{ii},jj,num_embeddings,remove_auto));
            distance_matrices(jj,kk,ii)=sum((conn_profilekk-conn_profileii).^2)/sqrt(sum(conn_profilekk.^2)*sum(conn_profileii.^2));
        end;
        
    end;
end;

function conn_profile=get_conn_profile(netmat_full,node_ind,num_embeddings,remove_auto)

from = (node_ind-1)*num_embeddings+1;
to = from+num_embeddings-1;
conn_profile=netmat_full(from:to,:);
if remove_auto
    conn_profile(:,from:to)=0;
end;
            
