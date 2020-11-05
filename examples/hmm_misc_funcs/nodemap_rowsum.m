function [ nodemap ] = nodemap_rowsum( state_netmats, S )

% [ nodemap ] = nodemap_rowsum( state_netmats, S )
%
% Computes sum of edge weights
% using state_netmats.state{k}.netmat_full
%
% state_netmats.state{k}.netmat is num_nodes x num_nodes x num_embeddings
% state_netmats.state{k}.netmat_full is (num_nodes x num_embeddings) x (num_nodes x num_embeddings)
%       where the (num_nodes x num_embeddings) dimension is in the order:
%       [node1,embedding1; node1,embedding2 ... node1,embeddingN;
%       [node2,embedding1; node2,embedding2 ... node2,embeddingN;
%       etc...
%
% nodemap{k} is num_nodes x 1 

try S.contrast=S.contrast; catch S.contrast='global_diff'; end;

num_embeddings=size(state_netmats.global.netmat,3);
num_nodes=size(state_netmats.global.netmat,1);

if size(state_netmats.global.netmat_full,1)==num_nodes,
    num_embeddings=1;
end;

remove_auto=S.remove_auto;

global_nodemap=compute_nodemap(state_netmats.global.netmat_full,num_embeddings,num_nodes,remove_auto);

NK=length(state_netmats.state);

for k=1:NK,
    nodemap(:,k)=compute_nodemap(state_netmats.state{k}.netmat_full,num_embeddings,num_nodes,remove_auto);
    switch S.contrast
        case 'none'
            % do nothing
        case 'global_diff'
            nodemap(:,k)=(nodemap(:,k)-global_nodemap')./global_nodemap'; 
    end;
end;

function nodemap=compute_nodemap(netmat,num_embeddings,num_nodes,remove_auto)

% collapse over node subdimension (e.g. freq or embedding)
nodemap=reshape(netmat,[num_embeddings,num_nodes,num_embeddings,num_nodes]);

nodemap=permute(sum(sum(abs(nodemap),1),3),[2 4 1 3]);

if remove_auto
    nodemap=nodemap-diag(diag(nodemap));
end;

nodemap=sum(nodemap,1);

