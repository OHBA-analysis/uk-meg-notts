function statePath = ABhmm_statepath(hmm)
% Decode HMM statepath
% statePath = ABhmm_statepath(hmm)
% AB 2013
  block = hmmdecode(hmm.data.Xtrain,size(hmm.data.Xtrain,1),hmm);
  statePath = abs(block(1).q_star);
end
