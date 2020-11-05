function [ nodemap res ] = nodemap_connprofile( state_netmats, S )

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

res=[];

num_embeddings=size(state_netmats.global.netmat,3);
num_nodes=size(state_netmats.global.netmat,1);
    
if size(state_netmats.global.netmat_full,1)==num_nodes,
    num_embeddings=1;
end;

NK=length(state_netmats.state);
        
clear nodemap;
for jj=1:num_nodes,
        
    for k=1:NK,
        % extract state profile for this node                 
        tmp=squash(full(get_conn_profile(state_netmats.state{k}.netmat_full,state_netmats.state{k}.netmat,jj,num_embeddings,S.mode)));                        
        tmp=tmp(tmp~=0);
        conn_profile(k,:)=tmp;
    end;
    
    mn=mean(conn_profile,1);
    for k=1:NK,   
        tmp=conn_profile(k,:);
        if 0
            thr=percentile(squash(conn_profile),75);            
            tmp(tmp<thr)=0;            
        end;
        
        nodemap(jj,k)=sqrt(sum((conn_profile(k,:)-mn).^2))./sqrt(sum(mn.^2));
       
        %nodemap(jj,k)=sum(tmp)/sum(mn);
    end;
end;

switch S.contrast_type
            
    case 'acope'
        nodemap=abs(nodemap);
    case 'cope'
        % do nothing
    otherwise
        error('unknown contrast_type');
        
end;

for k=1:NK,
    nodemap(:,k)=(nodemap(:,k)-mean(squash(nodemap(:,k))))/std(squash(nodemap(:,k)));
end;

function y=logit(x)

y=log(x./(1-x));

function conn_profile=get_conn_profile(netmat_full,netmat,node_ind,num_embeddings,mode)

from = (node_ind-1)*num_embeddings+1;
to = from+num_embeddings-1;

switch mode
    case 'all'
        conn_profile=netmat_full(from:to,:);
    case 'auto'
        conn_profile=netmat_full(from:to,from:to);
    case 'cross'
        conn_profile=netmat_full(from:to,:);
        conn_profile(:,from:to)=0;      
    otherwise
        error('Invalid mode');
end;
            
