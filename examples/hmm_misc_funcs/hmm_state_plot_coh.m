function hmm_state_plot_coh(state_netmats_mt, nodes, fc_type, state_for_subj_plots)


num_subj=length(state_netmats_mt);

cols={'r','g','b','m','c','y','r--','g--','b--','m--','c--','y--',};

if isfield(state_netmats_mt{1},'state')
    NK=length(state_netmats_mt{1}.state);
else
    NK=0;
end;

fs=state_netmats_mt{1}.global.spectramt.f;
clear legs;
legs=[];

abs_before_averaging=true;

% global
tot_subj_ntpts_global=0;
for ss=1:num_subj    
    try
        tot_subj_ntpts_global=tot_subj_ntpts_global+state_netmats_mt{ss}.global.spectramt.ntpts;
    catch
        state_netmats_mt{ss}.global.ntpts=1;
        tot_subj_ntpts_global=tot_subj_ntpts_global+state_netmats_mt{ss}.global.ntpts;
        
    end
end

clear psd_global coh_global;
for ss=1:num_subj
    psd_global(:,:,:,ss)=state_netmats_mt{ss}.global.spectramt.psd;
    coh_global(:,:,:,ss)=state_netmats_mt{ss}.global.spectramt.(fc_type);    
    if abs_before_averaging    % note: coherence is already real and postive
        psd_global(:,:,:,ss)=abs(psd_global(:,:,:,ss));
    end
    
    try
        psd_global(:,:,:,ss)=psd_global(:,:,:,ss)*sqrt(state_netmats_mt{ss}.global.spectramt.ntpts/tot_subj_ntpts_global);
        coh_global(:,:,:,ss)=coh_global(:,:,:,ss)*sqrt(state_netmats_mt{ss}.global.spectramt.ntpts/tot_subj_ntpts_global);
    catch
        psd_global(:,:,:,ss)=psd_global(:,:,:,ss)*sqrt(state_netmats_mt{ss}.global.ntpts/tot_subj_ntpts_global);
        coh_global(:,:,:,ss)=coh_global(:,:,:,ss)*sqrt(state_netmats_mt{ss}.global.ntpts/tot_subj_ntpts_global);
    end
end

tot_subj_ntpts=[];
for kk=1:NK
    tot_subj_ntpts(kk)=0;
    for ss=1:num_subj    
        try
            tot_subj_ntpts(kk)=tot_subj_ntpts(kk)+state_netmats_mt{ss}.state{kk}.spectramt.ntpts;
        catch
            state_netmats_mt{ss}.state{kk}.ntpts=1;

            tot_subj_ntpts(kk)=tot_subj_ntpts(kk)+state_netmats_mt{ss}.state{kk}.ntpts;
        end
    end
end

% statewise
figure;
for kk=1:NK
    
    legs{kk}=num2str(kk);
    clear psd coh;
    for ss=1:num_subj
        try
            psd(:,:,:,ss)=state_netmats_mt{ss}.state{kk}.spectramt.psd;
            coh(:,:,:,ss)=state_netmats_mt{ss}.state{kk}.spectramt.(fc_type);
            if abs_before_averaging  % note: coherence is already real and postive  
                psd(:,:,:,ss)=abs(psd(:,:,:,ss));
            end        

            if false
                % sanity check
                l=37;
                j=38;
                psdtmp=state_netmats_mt{ss}.state{kk}.spectramt.psd;

                % these next two should be the same
                cjl = abs(psdtmp(:,j,l)./sqrt(psdtmp(:,j,j) .* psdtmp(:,l,l)))
                state_netmats_mt{ss}.state{kk}.spectramt.coh(:,j,l)
            end

            try
                psd(:,:,:,ss)=psd(:,:,:,ss)*sqrt(state_netmats_mt{ss}.state{kk}.spectramt.ntpts/tot_subj_ntpts(kk));
                coh(:,:,:,ss)=coh(:,:,:,ss)*sqrt(state_netmats_mt{ss}.state{kk}.spectramt.ntpts/tot_subj_ntpts(kk));
            catch
                psd(:,:,:,ss)=psd(:,:,:,ss)*sqrt(state_netmats_mt{ss}.state{kk}.ntpts/tot_subj_ntpts(kk));
                coh(:,:,:,ss)=coh(:,:,:,ss)*sqrt(state_netmats_mt{ss}.state{kk}.ntpts/tot_subj_ntpts(kk));            
            end
        catch
            psd(:,:,:,ss)=nan;
            coh(:,:,:,ss)=nan; 
        end
    end

    try
        plot_cohs(abs(nanmean(psd,4)),abs(nanmean(coh,4)),nodes,fs,cols{kk},2,fc_type);
    catch
    end
end

plot_cohs(abs(mean(psd_global,4)),abs(mean(coh_global,4)),nodes,fs,'k',3,fc_type);       
legs{end+1}='g';

legend(legs);
subplot(2,2,2);title([mat2str(nodes)]);
    
if state_for_subj_plots>0
    
        k=state_for_subj_plots;
        figure;
        
        clear psd coh;
        for ss=1:num_subj
            psd(:,:,:,ss)=state_netmats_mt{ss}.state{k}.spectramt.psd;
            coh(:,:,:,ss)=state_netmats_mt{ss}.state{k}.spectramt.(fc_type);     
            ss2=rem(ss,length(cols))+1;
            plot_cohs(psd(:,:,:,ss),coh(:,:,:,ss),nodes,fs,cols{ss2},2,fc_type);
            hold on;
        end

        plot_cohs(nanmean(psd(:,:,:,1:num_subj),4),nanmean(coh(:,:,:,1:num_subj),4),nodes,fs,'k',3,fc_type);
        title(['State ' num2str(k)]);        
        
end

function legs=plot_cohs(psd, coh,nodes,fs,col,lw,fc_type)

if nargin<6
    lw=2;
end

if nargin<7
    fc_type='';
end

if size(psd,3)==1
    tmp=(psd(:,nodes(1)));
else
    tmp=(psd(:,nodes(1),nodes(1)));
end
subplot(2,2,1);plot(fs,tmp,col,'LineWidth',lw);ho;    
title(['psd[' mat2str(nodes(1)) ']']);

tmp=(coh(:,nodes(1),nodes(2)));
subplot(2,2,2);plot(fs,tmp,col,'LineWidth',lw);ho;    
title([fc_type  mat2str(nodes) ]);

tmp=(coh(:,nodes(2),nodes(1)));
subplot(2,2,3);plot(fs,tmp,col,'LineWidth',lw);ho;    
title([fc_type  mat2str(nodes) ]);

if size(psd,3)==1
    tmp=(psd(:,nodes(2)));
else
    tmp=(psd(:,nodes(2),nodes(2)));
end
subplot(2,2,4);plot(fs,tmp,col,'LineWidth',lw);ho;    
title(['psd[' mat2str(nodes(2)) ']']);

