function netstats = ABhmm_stats(statepath,fs,do_plots)
% Compute summary statistics from statepath. Specify sampling rate fs to 
% convert to time.
% netstats = ABhmm_stats(statepath,fs,do_plots)
% AB 2013

if isstruct(statepath) % probably an hmm structure so get the statepath
  if isfield(statepath,'fsample')
    fs = statepath.fsample;
  end
  block = hmmdecode(statepath.data.Xtrain,size(statepath.data.Xtrain,1),statepath);
  statepath = abs(block(1).q_star);
else
    statepath = statepath(:)';
  if nargin == 1
  fs = 1;
  end
end

for s = unique(statepath)
  

  lifetimes = diff(logical2epoch(statepath==s),[],2)./fs; lifetimes = lifetimes(lifetimes~=0);
  intervals = diff(logical2epoch(statepath~=s),[],2)./fs; intervals = intervals(intervals~=0);
  
  netstats(s).FractionalOccupancy =  sum(statepath==s) / length(statepath);
  netstats(s).nOccurrences = length(lifetimes);
  netstats(s).MeanIntervalLength = mean(intervals);
  netstats(s).MeanLifeTime = mean(lifetimes);
  netstats(s).LifeTimes = lifetimes;
  netstats(s).Intervals = intervals;  
  
  
end


if exist('do_plots','var') && do_plots
  
  fnames = fieldnames(netstats);
  fnames = fnames(1:4); %only plot mean stats
  for f = 1:length(fnames)
    ABsubplot(length(fnames),f), bar([netstats.(fnames{f})]), xlabel('State #'), title(fnames(f));
  end
  
  
  
end

end


function ev = logical2epoch(l,t)
% ev = logical2epoch(l,t)

if nargin < 2 
  t = 1:length(l);
end

onsets = t(diff([0 l])==  1);
offsets = t(diff([l 0])==-1);

ev = [onsets; offsets]';

%figure, plot(l), ho, scatter(onsets,ones(size(onsets)),'xg'), scatter(offsets,ones(size(offsets)),'xr')

end

