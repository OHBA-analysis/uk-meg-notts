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

if ~isfield(S,'demean')
    S.demean=1;
end
num_embeddings=size(state_netmats.global.netmat,3);
num_nodes=size(state_netmats.global.netmat,1);
    
if size(state_netmats.global.netmat_full,1)==num_nodes,
    num_embeddings=1;
end;

if isfield(state_netmats,'state')
    NK=length(state_netmats.state);
else
    NK=1;
    state_netmats.state{1}=state_netmats.global;
end;
        
res.conns=zeros(NK,num_nodes,num_nodes);
for k=1:NK,    
    for ii=1:num_nodes,            
        for jj=1:ii-1,
            conn=get_conn(state_netmats.state{k}.netmat_full,state_netmats.state{k}.netmat,ii,jj,num_embeddings);        
            
            % collapse over embeddings
            res.conns(k,ii,jj)=mean(abs(squash(conn)));
            res.conns(k,jj,ii)=res.conns(k,ii,jj);
        end;
    end;    
end;

clear nodemap;
for jj=1:num_nodes,
    % extract global profile for this node                 
    global_conn_profile=get_conn_profile(state_netmats.global.netmat_full,state_netmats.global.netmat,jj,num_embeddings,S.mode);
    if ~issparse(global_conn_profile)
        global_conn_profile=squeeze(mean(abs(global_conn_profile),1));
    else
        global_conn_profile=full(global_conn_profile);
        global_conn_profile=squash(global_conn_profile,global_conn_profile);
    end;     
    for k=1:NK,
        % extract state profile for this node                 
        conn_profile=get_conn_profile(state_netmats.state{k}.netmat_full,state_netmats.state{k}.netmat,jj,num_embeddings,S.mode);
        if ~issparse(conn_profile)
            conn_profile=squeeze(mean(abs(conn_profile),1));
        else
            conn_profile=full(conn_profile);
            conn_profile=squash(conn_profile,conn_profile);
        end;
        
        %res.connprofiles(jj,k,:)=conn_profile;
                
        switch S.contrast_type
            
            case 'coape'
                nodemap(jj,k)=median(conn_profile-global_conn_profile)./median(global_conn_profile);
            case 'coape'
                nodemap(jj,k)=median(conn_profile-global_conn_profile)./median(global_conn_profile);
            case 'mean' % areas where average connectivity is high
                nodemap(jj,k)=mean(conn_profile);
            case 'median' % areas where median connectivity is high
                nodemap(jj,k)=median(conn_profile);         
            case 'diffmean' % areas where overall connectivity increases most compared to global
                nodemap(jj,k)=mean(conn_profile)-mean(global_conn_profile);
            case 'diffmedian' % areas where overall connectivity increases most compared to global
                nodemap(jj,k)=(median(conn_profile)-median(global_conn_profile));
            case 'globalmedian'                
                nodemap(jj,k)=median(global_conn_profile);  
            case 'binding' % areas where variability in their edges is high over states (i.e. over time)
                nodemap_tmp(jj,k,:)=conn_profile;  
                nodemap(jj,k)=0;
            otherwise
                error('unknown contrast_type');
        end
    end
    
    % now perform any necessary operations over all states
    switch S.contrast_type
        case 'binding' % areas where variability in their edges is high over states (i.e. over time)
             nodemap(jj,:)=median(squeeze(std(nodemap_tmp(jj,:,:),[],2))); % std is over states, median is over conn profile (i.e. nodes/lags)
        otherwise
             % do nothing
    end
end

% normalise within state
for k=1:NK,
    if S.demean
        nodemap(:,k)=(nodemap(:,k)-mean(squash(nodemap(:,k))));
        nodemap(:,k)=(nodemap(:,k))/std(squash(nodemap(:,k)));
    end
end

if 1==0;
end

function y=logit(x)

y=log(x./(1-x));

function conn_profile=get_conn_profile(netmat_full,netmat,node_ind,num_embeddings,mode)

% returns num_embeddings x (num_embeddings x nnodes)

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
            
