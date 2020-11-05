function [W,omega,pred,residuals,cod] = MLmar (X,T,order,W,omega)

ndim = size(X,2);
N = length(T);

XX = []; Y = [];
for in=1:N
    if in==1, t0 = 0;
    else t0 = sum(T(1:in-1));
    end
    if order>0
        XX0 = zeros(T(in)-order,order*ndim);
        for i=1:order
            XX0(:,(1:ndim) + (i-1)*ndim) = X(t0+order-i+1:t0+T(in)-i,:);
        end;
        XX = [XX; XX0];
    end
    Y = [Y; X(t0+order+1:t0+T(in),:)];
end


if nargin<4 || isempty(W)
    W = zeros(order*ndim,ndim);
    
    if order>0
        for n=1:ndim
            W(:,n) = (XX' * XX) \ XX' * Y(:,n);
        end
        pred = XX * W;
    else
        pred = zeros(sum(T),ndim);
    end
    
else
    pred = XX * W;
    
end


residuals = Y - pred;
SSE = sum ( residuals.^2 );
omega = (residuals' * residuals ) / (size(Y,1)-1) ;
cod = 1 - SSE ./ sum( (Y - repmat( mean(Y),size(Y,1),1) ).^2 );


