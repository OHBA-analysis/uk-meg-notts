function KL = ABhmm_kl(hmm)
% Compute KL divergence of HMM states
% AB 2013

KL = zeros(hmm.K);

for i = 1:hmm.K
  for j = 1:hmm.K
    
    ci = hmm.state(i).Cov;
    cj = hmm.state(j).Cov;
    
    KL(i,j) = trace(pinv(cj)*ci)-log(det(ci)/det(cj))-size(cj,1);
    
  end
end

KL(logical(eye(hmm.K))) = nan;

end
