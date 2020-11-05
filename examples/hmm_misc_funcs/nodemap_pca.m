function [ nodemap, res ] = nodemap_pca( state_netmats, S )

% [ nodemap ] = nodemap_pca( state_netmats, S )
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
        
clear nodemap res;

tmp_global=reshape(state_netmats.global.netmat,[num_nodes,num_nodes*num_embeddings]);

%tmp_global=normalise(tmp_global,2);
%pcadim = 1;
%[allsvd,Mglobal] = eigdec(cov(tmp_global'),pcadim);

for kk=1:NK,
    tmp=reshape(state_netmats.state{kk}.netmat,[num_nodes,num_nodes*num_embeddings]);
    tmp=tmp-tmp_global;
    
    tmp=normalise(tmp,2);
    pcadim = 1;
    [allsvd,M] = eigdec(cov(tmp'),pcadim);
    nodemap(:,kk) = abs(M);
    tmp2=M' * tmp;
    res.embed_map(:,:,kk) = reshape(tmp2,[num_nodes,num_embeddings]); % should be 1 x ntpts
end;


for k=1:NK,
    %nodemap(:,k)=(nodemap(:,k)-mean(squash(nodemap(:,k))))/std(squash(nodemap(:,k)));
    %nodemap(:,k)=(nodemap(:,k))/std(squash(nodemap(:,k)));
end;