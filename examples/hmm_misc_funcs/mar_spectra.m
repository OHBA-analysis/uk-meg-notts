function [fit,freqs] = mar_spectra (A,covm,Hz,Nf,minHz,maxHz)
% Get spectral estimates from MAR model
%
% A: MAR autoregression coefficients
% covm: covariance matrix of the residuals
% ns:    samples per second
% Nf:    number of frequencies to plot
% minHz - maxHz: band of frequencies to look at 
%
% The returned structure will have the following fields specified:
%
% .PSD     [Nf x d x d] Power Spectral Density matrix
% .iPSD     [Nf x d x d] Inverse Power Spectral Density matrix
% .C     [Nf x d x d] Coherence matrix
% .pC    [Nf x d x d] Partial Coherence matrix
% .dtf   [Nf x d x d] Kaminski's Directed Transfer Function matrix
% .pdc   [Nf x d x d] Baccala's Partial Directed Coherence
% .L     [Nf x d x d] Phase matrix
% .f     [Nf x 1] Frequency vector
% .ns    Sample rate
%
% dtf(f,i,j) is the DTF at frequency f from signal j to signal i
% pdc(f,i,j) is the PDC at frequency f from signal j to signal i
%

if nargin<3
    Nf = 512;
end

order=size(A,1) / size(A,2);
ndim=size(A,2);

freqs = (0:Nf-1)*( (maxHz - minHz) / (Nf-1)) + minHz;
w=2*pi*freqs/Hz;

W = zeros(order,ndim,ndim);
for i=1:order,
    W(i,:,:) = A((1:ndim) + (i-1)*ndim,:);
end
prec = inv(covm);
fit = {};

% Get Power Spectral Density matrix and DTF for state K
for ff=1:Nf,
    af_tmp=eye(ndim);
    for o=1:order,
        af_tmp=af_tmp - squeeze(W(o,:,:)) * exp(-1i*w(ff)*o);
    end
    iaf_tmp=inv(af_tmp);
    PowSpD(ff,:,:) = iaf_tmp * covm * iaf_tmp';
    iPowSpD(ff,:,:) = inv(permute(PowSpD(ff,:,:),[3 2 1]));
    fit.PSD(ff,:,:) = real(PowSpD(ff,:,:)).^2;
    fit.iPSD(ff,:,:) = real(iPowSpD(ff,:,:)).^2;
    
    % Get DTF and PDC
    for nn=1:ndim,
        prec_nn=1/sqrt(covm(nn,nn));
        for n=1:ndim,
            % DTF uses iaf_tmp and normalises wrt rows (to=sink)
            fit.dtf(ff,nn,n)=abs(iaf_tmp(nn,n))/sqrt(iaf_tmp(nn,:)*iaf_tmp(nn,:)');
            % PDC uses af_tmp and normalises wrt columns (from=source)
            fit.pdc(ff,nn,n) = prec_nn * abs(af_tmp(nn,n))/sqrt(abs(af_tmp(:,n)'*prec*af_tmp(:,n)));
        end
    end
    
end

% Get Coherence and Phase
for nn=1:ndim,
    for n=1:ndim,
        rkj=PowSpD(:,nn,n)./(sqrt(PowSpD(:,nn,nn)).*sqrt(PowSpD(:,n,n)));
        %fit.C(:,nn,n)=abs(rkj);
        fit.C(:,nn,n)=rkj;
        fit.pC(:,nn,n)=-iPowSpD(:,nn,n)./(sqrt(iPowSpD(:,nn,nn)).*sqrt(iPowSpD(:,n,n)));
        l=atan(imag(rkj)./real(rkj));
        fit.L(:,nn,n)=atan(imag(rkj)./real(rkj));
    end
end

end





