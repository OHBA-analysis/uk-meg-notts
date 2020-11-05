function [ nodemap ] = nodemap_connprofile( state_netmats, S )

% [ nodemap ] = nodemap_connprofile( state_netmats, S )
%
% Computes connectivity profile distance from average over all states
% using state_netmats.state{k}.netmat_full
%
% state_netmats.state{k}.netmat_full is (num_nodes x num_embeddings) x (num_nodes x num_embeddings)
%       where the (num_nodes x num_embeddings) dimension is in the order:
%       [node1,embedding1; node1,embedding2 ... node1,embeddingN;
%       [node2,embedding1; node2,embedding2 ... node2,embeddingN;
%       etc...
%
% nodemap{k} is num_nodes x 1 

num_embeddings=size(state_netmats.global.netmat,3);
num_nodes=size(state_netmats.global.netmat,1);
    
if size(state_netmats.global.netmat_full,1)==num_nodes,
    num_embeddings=1;
end;

NK=length(state_netmats.state);
        
clear nodemap;
for jj=1:num_nodes,
    % extract global profile for this node                 
    global_conn_profile=squash(get_conn_profile(state_netmats.global.netmat_full,jj,num_embeddings,S.mode));
        
    for k=1:NK,
        % extract state profile for this node                 
        conn_profile=squash(get_conn_profile(state_netmats.state{k}.netmat_full,jj,num_embeddings,S.mode));

        nodemap(jj,k)=sqrt(sum((conn_profile-global_conn_profile).^2))./sqrt(sum((global_conn_profile).^2));
    end;
end;

function conn_profile=get_conn_profile(netmat_full,node_ind,num_embeddings,mode)

from = (node_ind-1)*num_embeddings+1;
to = from+num_embeddings-1;

switch mode
    case 'all'
        conn_profile=netmat_full(from:to,:);
        conn_profile(:,from:to)=0;
    case 'auto'
        conn_profile=netmat_full(from:to,from:to);
    case 'cross'
        conn_profile=netmat_full(from:to,:);
        conn_profile(:,from:to)=0;        
end;
            
