function stat = ABhmm_regress(hmm,voxeldata,use_abs,mode)
% Fits the HMM statepath as a regressor on the voxelwise data  
% stat = ABhmm_regress(hmm,voxeldata,use_abs,mode)
% mode can be 'cope','tstat','corr','pcorr'
% AB 2013

if ~exist('use_abs','var') || isempty(use_abs)
  use_abs = 0;
end

if ~exist('mode','var')
  mode = 'cope';
end




if isstruct(hmm)
  block = hmmdecode(hmm.data.Xtrain,size(hmm.data.Xtrain,1),hmm);
  statePath = abs(block(1).q_star);
  Nstates = hmm.K;
  if size(voxeldata,2) ~= size(hmm.data.Xtrain,1)
    voxeldata = voxeldata';
  end
  if isfield(hmm,'whiteningMatrix')
    Nvoxels = size(hmm.whiteningMatrix,2);
  else
    Nvoxels = size(voxeldata,1);
  end
else
  statePath = hmm;
  Nstates = length(unique(statePath));
  if size(voxeldata,2) ~= length(statePath)
    voxeldata = voxeldata';
  end
  Nvoxels = size(voxeldata,1);
end




% Regress Viterbi path onto wholebrain results
for i=1:Nstates
    con{i}=(1:Nstates==i)';
    %con{i}=(1:Nstates==i)' - (1:Nstates~=i)'/(Nstates-1);
end


cope    = zeros(Nvoxels,Nstates);
varcope = zeros(Nvoxels,Nstates);
c       = zeros(Nvoxels,Nstates); 

x = zeros(length(statePath),Nstates);
for s = 1:Nstates
  x(:,s) = double(statePath == s);
end

if strcmp(mode,'pcorr')
  x = devar(x,1);
else
  x = demean(x,1);
end

pinvxtx = pinv(x'*x);
pinvx = pinv(x);
  


for v=1:Nvoxels
  
  if isempty(voxeldata);
    iwm = pinv(hmm.whiteningMatrix);
    voxeldata_v = iwm(v,:)*hmm.data.Xtrain';
    %voxeldata_v = pinv(hmm.whiteningMatrix(:,v))*hmm.data.Xtrain';
  else
    voxeldata_v = voxeldata(v,:);
  end
  
  if use_abs
    if strcmp(mode,'pcorr')
      y = normalise(abs(hilbert(voxeldata_v))');
    else
      y = demean(abs(hilbert(voxeldata_v))');
    end
  else
    if strcmp(mode,'pcorr')
      y = normalise(voxeldata_v)';
    else
      y = demean(voxeldata_v)';
    end
  end
    
  if strcmp(mode,'corr')
    c(v,:) = corr(x,repmat(y,1,Nstates));
  else
    [cope(v,:) varcope(v,:)] = glm_fast_for_meg(y,x,pinvxtx,pinvx,con,0);
  end

end

tstat = cope ./ sqrt(varcope);

switch mode
  case {'cope','pcorr'}
    stat = cope;
  case 'tstat'
    stat = tstat;
  case 'corr'
    stat = c;
end
