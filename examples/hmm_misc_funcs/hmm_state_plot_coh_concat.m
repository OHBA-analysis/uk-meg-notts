function legs=hmm_state_plot_coh_concat(state_netmats_mt, nodes, fc_type, do_diff, print_plot)

% do_diff=1 

cols=get_cols();

if isfield(state_netmats_mt{1},'state')
    NK=length(state_netmats_mt{1}.state);
else
    NK=0;
end;

if nargin<5
    print_plot=0;
end

fs=state_netmats_mt{1}.global.spectramt.f;
clear legs;
legs=[];

% global
clear psd_global coh_global;
ss=1;

psd_global(:,:,:,ss)=abs(state_netmats_mt{ss}.global.spectramt.psd);
coh_global(:,:,:,ss)=abs(state_netmats_mt{ss}.global.spectramt.(fc_type));

% statewise
figure;
for k=1:NK
    
    legs{k}=num2str(k);
    psd(:,:,:,ss)=abs(state_netmats_mt{ss}.state{k}.spectramt.psd);
    coh(:,:,:,ss)=abs(state_netmats_mt{ss}.state{k}.spectramt.(fc_type));

    if do_diff
        plot_cohs(psd-psd_global,coh-coh_global,nodes,fs,cols{k},2,fc_type,print_plot);       
    else
        plot_cohs(psd,coh,nodes,fs,cols{k},2,fc_type,print_plot);
    end;
    
end;

if ~do_diff
    plot_cohs(psd_global,coh_global,nodes,fs,'k',3,fc_type,print_plot);       
    legs{end+1}='g';
end;

if ~print_plot
    legend(legs);    
end

if print_plot
    subplot(2,2,1);plot4paper('freqs (Hz)','');
    %set(gca,'YtickLabels',[]);
    subplot(2,2,2);plot4paper('freqs (Hz)','')
    %set(gca,'YtickLabels',[]);
    subplot(2,2,3);plot4paper('freqs (Hz)','')
    %set(gca,'YtickLabels',[]);
    subplot(2,2,4);plot4paper('freqs (Hz)','')
    %set(gca,'YtickLabels',[]);
    
end

function legs=plot_cohs(psd, coh,nodes,fs,col, lw, fc_type,print_plot)

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
if ~print_plot,title(['psd[' mat2str(nodes(1)) ']']);end

tmp=(coh(:,nodes(1),nodes(2)));
subplot(2,2,2);plot(fs,tmp,col,'LineWidth',lw);ho;    
if ~print_plot, title([fc_type  mat2str(nodes) ]); end

tmp=(coh(:,nodes(2),nodes(1)));
subplot(2,2,3);plot(fs,tmp,col,'LineWidth',lw);ho;    
if ~print_plot,title([fc_type  mat2str(nodes) ]);end

if size(psd,3)==1
    tmp=(psd(:,nodes(2)));
else
    tmp=(psd(:,nodes(2),nodes(2)));
end
subplot(2,2,4);plot(fs,tmp,col,'LineWidth',lw);ho;    
if ~print_plot,title(['psd[' mat2str(nodes(2)) ']']);end

