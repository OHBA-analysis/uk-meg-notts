function tstats_fnames=event_related_average(S)

D=S.D;
baseline=S.baseline;
results_fname=S.results_fname;

for condnum = 1:length(D.condlist),

    clear env;

    trials = D.indtrial(D.condlist{condnum},'good'); 

    for trl=1:length(trials),           

        tbad = all(badsamples(D,':',':',trials(trl)));
        samples2use = find(~tbad);

        env(:,samples2use,trl) = D(:,samples2use,trials(trl)); %#ok - SPM doesn't like logical indexing
    end;

    ers=mean(env,3);
    stders=std(env,[],3).^2;

    % baseline correct
    mn=mean(ers(:,find(D.time(samples2use)>-8 & D.time(samples2use)<-4)),2);
    ers(:,samples2use)=ers(:,samples2use)-repmat(mn,1,length(samples2use));

    ts=ers./sqrt(stders/size(env,3));

    gridstep=getmasksize(size(ts,1));
    tstats_fnames{condnum}=nii_quicksave(ts,[results_fname '_tstat_' D.condlist{condnum}],gridstep,gridstep);

end;

