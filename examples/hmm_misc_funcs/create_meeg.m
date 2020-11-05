function D =create_meeg(data_new, fname, Fs)

nregions=size(data_new,2);

chan_labels=[];
for reg = 1:nregions
    chan_labels{reg}=['chan' num2str(reg)];        
end

ftdata=[];
ftdata.label=chan_labels; % the channel labels
tbins = 1:size(data_new,1);
ftdata.time{1} = tbins/Fs; % the timeaxis [1*Ntime double] per trial
ftdata.trial{1} = data_new'; % the numeric data [Nchans*Ntime] per trial

D=spm_eeg_ft2spm(ftdata, fname);
D.save;

end

