function fit = fullcov(data,Gamma,options)

fit=[];

K=size(Gamma,2);

[aa bb]=max(Gamma');
    
for kk = 1:K

    disp(['Computing for state ' num2str(kk)]);

    inds = logical(bb==kk)';

    data_in=data(inds)';
            
    [fit.state{kk}]=netmat_cov(data_in,options);

end