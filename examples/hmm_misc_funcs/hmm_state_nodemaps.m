function [ res ] = hmm_state_nodemaps( S )

% [ res ] = hmm_state_nodemaps( S )

try S.output_separate_embeddings=S.output_separate_embeddings; catch S.output_separate_embeddings=0; end;
try mask_fname=S.mask_fname; catch mask_fname=[]; end
try S.netmat_node_method_options=S.netmat_node_method_options; catch S.netmat_node_method_options=[]; end
%try S.subj_ind=S.subj_ind; catch S.subj_ind=1; end

nsubs=length(S.state_netmats);

clear res_netmat_node_method;
for subnum = 1:nsubs
    subnum
    % call passed in method
    [nodemaps_tmp res_netmat_node_method{subnum}]=feval(S.netmat_node_method,S.state_netmats{subnum},S.netmat_node_method_options);
    
    if subnum==1
        nodemaps=nodemaps_tmp;
    else
        nodemaps=nodemaps+nodemaps_tmp;
    end;
end;
nodemaps=nodemaps/nsubs;

res.res_netmat_node_method=res_netmat_node_method;

% save nii file
S2=[];
S2.interp='nearestneighbour';
S2.mask_fname=mask_fname;
S2.output_spat_res=2; %mm

node_maps_fname=S.node_maps_fname;
% 
% if isfield(S.state_netmats{1}.parcellation,'S') && isfield(S.state_netmats{1}.parcellation.S,'hcp_sourcemodel3d') 
%     res.fname=hcp_nii_parcel_quicksave(nodemaps,S.state_netmats{1}.parcellation.assignments,node_maps_fname,S.state_netmats{1}.parcellation.S.hcp_sourcemodel3d, S2.output_spat_res, S2);
% else
%     res.fname=ROInets.nii_parcel_quicksave(nodemaps,S.state_netmats{1}.parcellation.assignments,node_maps_fname,S2);
% end

map = parcellation2map(nodemaps,S.parcellation.file,S.parcellation.mask);
writenii(map,node_maps_fname,S.parcellation.mask);
res.fname=node_maps_fname;
     
%res.fname=ROInets.nii_parcel_quicksave(nodemaps,S.state_netmats{1}.parcelAssignments,node_maps_fname,S2);
    
res.nodemaps=nodemaps;

if isfield(S.state_netmats{1},'state')
    NK=length(S.state_netmats{1}.state);
else
    NK=0;
end;
                                                                                    
if S.output_separate_embeddings
    
    clear nodemaps_byfreq;
    for ff=1:size(S.state_netmats.state{1}.netmat,3),

        for subnum = 1:nsubs

            state_netmats=S.state_netmats{subnum};

            state_netmats.global.netmat=state_netmats.global.netmat(:,:,ff);
            state_netmats.global.netmat_full=S.state_netmats{subj_ind}.global.netmat;

            % setup up state_netmats for freq  
            for k=1:NK, 
                state_netmats.state{k}.netmat=state_netmats.state{k}.netmat(:,:,ff);
                state_netmats.state{k}.netmat_full=S.state_netmats{subj_ind}.state{k}.netmat;
            end;

            S.netmat_node_method_options.output_separate_embeddings=S.output_separate_embeddings;
            
            nodemaps_byfreq(:,:,ff)=nodemaps_byfreq(:,:,ff)+feval(S.netmat_node_method,state_netmats,S.netmat_node_method_options);
            
        end;

        nodemaps_byfreq(:,:,ff)=nodemaps_byfreq(:,:,ff)/nsubs;
    end;
    
    for k=1:NK, 
                        
        S2=[];
        S2.interp='nearestneighbour';
        S2.mask_fname=mask_fname;
        node_maps_fname=[S.node_maps_fname '_k' num2str(k)];        
        res.fname_byfreq{k}=ROInets.nii_parcel_quicksave(permute(nodemaps_byfreq(:,k,:),[1 3 2]),S.state_netmats.parcelAssignments,node_maps_fname,S2);

    end
    
    res.fs=S.state_netmats.global.spectramt.f;
end

end

