function ABhmm_plotFO(statepath,fs,win,mode,t0,bc)
% Plots the HMM statepath
% ABhmm_plotFO(statepath,fs,win)
% AB 2013

%%% MWW added
if nargin < 4
  mode = 'separated'; % stack all states on separate lines on y-axis
end
if nargin < 5
  t0 = 0; % start time point in secs
end
if nargin < 6
  bc = []; % baseline correct epoched states using this time range
end
%%%

k = max(statepath);

clf, hold on
temp_fig = figure('visible','off'); col = colormap(temp_fig,'lines'); close(temp_fig)

leg=[];
for s = 1:k
  %plot((1:length(statepath))/fs, 0.8*double(statepath==s)+ s - 0.4,'color',col(s,:));
  
  %% MWW to handle epoch data
  tmp=double(logical(statepath==s));
  if size(tmp,1)>1,
      %epoched data
      tmp=mean(tmp,2);
  end;
  
  FO = conv(tmp,rectwin(fs*win),'same')./ (fs*win);
  ts=linspace(t0,t0+length(FO)/fs,length(FO));

  %% MWW 
  if ~isempty(bc)
    FO=FO-mean(FO(intersect(ts>bc(1),ts < bc(2))))
  end;
  
  %% MWW to handle epoch data  
  if strcmp(mode,'separated'),      
    plot(ts,0.8*FO + s - 0.4,'color',col(s,:));
    xlabel('Time (s)'); ylabel('State #')
  else
    plot(ts,FO,'color',col(s,:));  
    xlabel('Time (s)'); ylabel('FO')
    leg=[leg; 's' num2str(s)];
    
  end;
end

%% MWW to handle epoch data  
if ~strcmp(mode,'separated'),
    legend(leg);
end;

hold off

end
