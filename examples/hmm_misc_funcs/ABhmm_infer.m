function hmm = ABhmm_infer(hmmdata,nstates,nreps,varargin)
% Runs an HMM with nstates multiple times
% hmm = ABhmm_infer(hmmdata,nstates,nreps,varargin)
% varargin: {'constrain_mean','do_stats'}
% AB 2013

FrEn = Inf;
hmm = struct;
hmm.FrEn  = [];
hmm.stats = [];

if ~all(cellfun(@isstr,varargin))
  error('varargin must be a cell array of strings')
end

if any(ismember(varargin,'constrain_mean'))
  constrain_mean = 1;
else
  constrain_mean = 0;
end

if any(ismember(varargin,'do_stats'))
  do_stats = 1;
  hmm.stats = [];
else
  do_stats = 0;
end



for n = 1:nreps
  
 % try % sometimes it fails...
    hmm_new                = struct('K',nstates);
    hmm_new                = hmminit(hmmdata,hmm_new,'full');
    hmm_new.train.cyc      = 200;
    hmm_new.obsmodel       = 'Gauss';
    hmm_new.train.obsupdate= ones(1,hmm_new.K);
    hmm_new.train.init     = 1;

    if constrain_mean
      disp('Constraining mean')
      for k=1:length(hmm_new.state)
        hmm_new.state(k).Mu = 0*hmm_new.state(k).Norm_Mu;
        hmm_new.state(k).constraint = 1;
      end
    else
      disp('Not constraining mean')
    end
    
    [hmm_new,FrEn_new] = hmmtrain(hmmdata,length(hmmdata),hmm_new);
    
    if FrEn_new < FrEn
      FrEn = FrEn_new;
      % copy persistent info to new hmm struct:
      hmm_new.FrEn = hmm.FrEn;
      hmm_new.stats = hmm.stats;
      % replace hmm struct with new
      hmm = hmm_new;
    end
    hmm.FrEn = [hmm.FrEn FrEn_new];
    
    hmm.statepath   = ABhmm_statepath(hmm_new); 
    
    if do_stats % save stats from this run
      hmmstats.state       = hmm_new.state; 
      hmmstats.transitions = hmm_new.P;
      hmmstats.statepath   = hmm.statepath;
      hmmstats.KL          = ABhmm_kl(hmm_new);
      hmm.stats = [hmm.stats hmmstats];
    end

    
 % end
  disp(['Elapsed time = ' num2str(cputime)])
end


end
