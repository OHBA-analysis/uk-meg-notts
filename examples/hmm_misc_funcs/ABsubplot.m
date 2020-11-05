function [i j] = ABsubplot(N,n)
i = ceil(sqrt(N));
j = ceil(N/i);

if nargin > 1
  subplot(i,j,n)
end

end
