function FO = ABhmm_smoothedFO(hmm,winsize,resamp,disp)
% ABhmm.statepath_smoothedFO(hmm,winsize)
% winsize in samples
% resamp 1 or 0 to resample or not
% disp 1 or 0 to plot


FO = zeros(length(unique(hmm.statepath)),length(osl_movavg(double(hmm.statepath==1),[],winsize,0.75,resamp)));

list=unique(hmm.statepath);
for i = 1:size(FO,1)
  
  FO(i,:) = osl_movavg(double(hmm.statepath==i),[],winsize,0.75,resamp);
end
  

if disp
  col = colormap(gcf,'lines'); ho
  for i = 1:size(FO,1)
    plot(1:size(FO,2), 0.8*FO(i,:)./max(FO(i,:)) + i - 0.4,'color',col(i,:));
  end  
end
  
end