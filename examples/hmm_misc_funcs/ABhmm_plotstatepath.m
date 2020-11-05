function ABhmm_plotstatepath(statepath,fs,subj_inds,ax)
% Plots the HMM statepath
% ABhmm_plotstatepath(statepath,fs)
% AB 2013

if nargin < 2
  fs = 1;
end

if nargin < 4
  ax = gca;
end

hold(ax,'on')

k = max(statepath);

col = colormap(ax,'lines');

for s = 1:k
  plot((1:length(statepath))/fs, 0.8*double(statepath==s)+ s - 0.4,'color',col(s,:));
end

if nargin > 2
  yl = get(gca,'ylim');
  stem(find(diff(subj_inds))/fs,yl(2)*ones(size(find(diff(subj_inds)))),'--k','marker','none');
  set(gca,'ylim',yl);
end

if fs ~= 1, xlabel('Time (s)'); else xlabel('samples'); end
ylabel('State #')

end
