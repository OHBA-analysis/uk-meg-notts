function plot_hmm_state_netmats( state_netmats, embed_ind, clim, subj_ind)

%  plot_hmm_state_netmats( state_netmats, embed_ind, clim, subj_ind)
% 
% state_netmats.state{k}.netmat is num_nodes x num_nodes x num_embeddings
% state_netmats.state{k}.netmat_full is (num_nodes x num_embeddings) x (num_nodes x num_embeddings)
% embed_ind should be emdedding dimension index, or 'full' to show
%       netmat_full (default)

if nargin<2
    embed_ind='full'; % plot full
end;
if nargin<3
    clim=[]; % plot full
end;
if nargin<4
    subj_ind=1;
end;

if isfield(state_netmats{1},'state')
    NK=length(state_netmats{1}.state);
else
    NK=0;
end;

figure;
for k = 1:NK
        
    switch embed_ind
        case 'full'
            mat=state_netmats{subj_ind}.state{k}.netmat_full;
        case 'mean_abs'
            mat=mean(abs(state_netmats{subj_ind}.state{k}.netmat),3);
        case 'mean'
            mat=mean((state_netmats{subj_ind}.state{k}.netmat),3);        
        otherwise
            embed_ind_num=str2num(embed_ind);
            mat=state_netmats{subj_ind}.state{k}.netmat(:,:,embed_ind_num);    
    end;
    
    mat=mat-diag(diag(mat));

    subplot(3,ceil((NK+1)/3),k);
    if isempty(clim)
        imagesc(mat);
    else
        imagesc(mat,clim);
    end;
    colorbar
end;

switch embed_ind
    case 'full'
        mat=state_netmats{subj_ind}.global.netmat_full;
    case 'mean_abs'
        mat=mean(abs(state_netmats{subj_ind}.global.netmat),3);
    case 'mean'
        mat=mean((state_netmats{subj_ind}.global.netmat),3);        
    otherwise
        embed_ind_num=str2num(embed_ind);
        mat=state_netmats{subj_ind}.global.netmat(:,:,embed_ind_num);    
end;

mat=mat-diag(diag(mat));

subplot(3,ceil((max(1,NK)+1)/3),NK+1);
if isempty(clim)
    imagesc(mat);
else
    imagesc(mat,clim);
end;
colorbar


end

