function hmm_state_plot_coh_new(state_netmats_mt, nodes, fc_type, do_diff, state_for_subj_plots)

% do_diff=1 does within subj paired t-test style comparisons with global PSDs

cols={'r','g','b','m','c','b--','k--','g--','r--'};

if isfield(state_netmats_mt{1},'state')
    NK=length(state_netmats_mt{1}.state);
else
    NK=0;
end;

fs=state_netmats_mt{1}.global.spectramt.f;
clear legs;
legs=[];

nsubs=length(state_netmats_mt);

% global
clear psd_global coh_global;
for ss=1:nsubs
    psd_global(:,:,:,ss)=abs(state_netmats_mt{ss}.global.spectramt.psd);
    coh_global(:,:,:,ss)=abs(state_netmats_mt{ss}.global.spectramt.(fc_type));
end;

% statewise
clear psd coh;
for k=1:NK
    
    legs{k}=num2str(k);
    for ss=1:nsubs
        psd(:,:,:,ss,k)=abs((state_netmats_mt{ss}.state{k}.spectramt.psd));
        coh(:,:,:,ss,k)=abs(state_netmats_mt{ss}.state{k}.spectramt.(fc_type));
    end;
   
end;

figure;
for k=1:NK
    inds=setdiff(1:NK,k);
    for ss=1:nsubs
        psd_global_minusk=mean(psd(:,:,:,ss,inds),5);
        coh_global_minusk=mean(coh(:,:,:,ss,inds),5);
        psd_diff(:,:,:,ss,k)=psd(:,:,:,ss,k)-psd_global_minusk;
        coh_diff(:,:,:,ss,k)=coh(:,:,:,ss,k)-coh_global_minusk;
        
    end;
    
    if do_diff
        plot_cohs(mean(psd_diff(:,:,:,:,k),4),mean(coh_diff(:,:,:,:,k),4),nodes,fs,cols{k});       
    else
        plot_cohs(mean((psd(:,:,:,:,k)),4),mean((coh(:,:,:,:,k)),4),nodes,fs,cols{k});
    end;
end;

if ~do_diff
    plot_cohs(mean((psd_global),4),mean((coh_global),4),nodes,fs,'k',2);    
    legs{end+1}='g';
end;

legend(legs);
subplot(2,2,2);title([mat2str(nodes)]);

if state_for_subj_plots>0
    
        k=state_for_subj_plots;
        figure;
        for ss=1:nsubs
            plot_cohs(psd(:,:,:,ss,k),coh(:,:,:,ss,k),nodes,fs,cols{ss});
            hold on;
        end;

        plot_cohs(mean(psd(:,:,:,:,k),4),mean(coh(:,:,:,:,k),4),nodes,fs,'k',2);
        title(['State ' num2str(k)]);
end;

function legs=plot_cohs(psd, coh,nodes,fs,col, lw)

if nargin<6
    lw=1;
end;

tmp=(psd(:,nodes(1),nodes(1)));
subplot(2,2,1);plot(fs,tmp,col,'LineWidth',lw);ho;    

tmp=(coh(:,nodes(1),nodes(2)));
subplot(2,2,2);plot(fs,tmp,col,'LineWidth',lw);ho;    

tmp=(coh(:,nodes(2),nodes(1)));
subplot(2,2,3);plot(fs,tmp,col,'LineWidth',lw);ho;    

tmp=(psd(:,nodes(2),nodes(2)));
subplot(2,2,4);plot(fs,tmp,col,'LineWidth',lw);ho;    


