function fit = spectramt(data,Gamma,options)

fit=[];

options.type='coh';
options.full_type='full';
options.var_normalise=true;
options.reg=.5;

K=size(Gamma,2);

[aa bb]=max(Gamma');

for kk = 1:K

    disp(['Computing for state ' num2str(kk)]);

    inds = logical(bb==kk)';

    data_in=data(inds);
            
    [fit.state{kk}]=netmat_spectramt(data_in,options);

end