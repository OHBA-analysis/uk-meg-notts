function [PLV, PSD] = PLI_band2(data, Fs, windowWidth, freqBands, doLeakageCorrection,doVarianceNormalisation)
%PLI  weighted phase lag index for pre-filtered data
%     use sym orthogonalisation within windows
%
%	computes:
%     Phase Lag Index
%     weighted Phase Lag Index
%     Phase Locking Value
%     Coherence
%     Imaginary Coherence
%     Phase Slope Index
%     Mutual Information of the phase
%
%    On continuous data, at sampling frequecny Fs, with window Width WS
%    seconds, looking only in frequency band FREQBAND for spectral
%    measures, and potentially applying leakage correction and z-transform
%    with data surrogates. 


%	References:
%	
%	

%	Copyright 2015 OHBA
%	This program is free software: you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%	
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%	
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.


%	$LastChangedBy$
%	$Revision$
%	$LastChangedDate$
%	Contact: giles.colclough@magd.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 09-May-2015 11:11:51

if nargin < 3 || ~exist('windowWidth', 'var') || isempty(windowWidth),
    windowWidth = 5; %s
end

%% Setup
[nSensors, nSamples] = size(data);

% we split resting data into windows, in a similar manner to Welch's method, to do averages.
windInd        = get_windows(nSamples, windowWidth, Fs);
windSampleSize = windowWidth * Fs;
nfft           = (windSampleSize) + mod(windSampleSize,2); % make even
nWindows       = ROInets.cols(windInd);
nFreqs=length(freqBands);

if nWindows>0
    % Declare memory
    PLVw(nSensors, nSensors, nWindows, nFreqs)          = 0;          
    ppsd=[];
    cpsd=[];

    %% Run calculation loop
    for iW = 1:nWindows,
        fprintf('Calculating PLI for window %d out of %d. \n', iW, nWindows);
        % cut data from window
        dataToUse  = ROInets.demean(data(:,windInd(:,iW)).',1);

        for iF=1:nFreqs,
            instabilityfix = 'reduce';
            dir='twopass';
            type = 'but';
            N = 5;
            Fp=freqBands{iF};
            dataToUse_filt = ft_preproc_bandpassfilter(dataToUse',Fs,Fp,N,type,dir,instabilityfix);
            dataToUse_filt=dataToUse_filt';

            % orthogonalise and variance normalise
            if doLeakageCorrection,
                dataToUse_filt  = ROInets.remove_source_leakage(dataToUse_filt,  'symmetric').';
            end%if

            if doVarianceNormalisation
                dataToUse_filt  = GC_normalise_vectors(dataToUse_filt, 1);
            end;

            % run pair-wise metrics
            for iSensor = 1:nSensors,

                iData = dataToUse_filt(:,iSensor);
                if iSensor>1,
                    thetaX=[];
                    for jSensor = 1:iSensor-1,

                            [PLVw(iSensor, jSensor, iW, iF), ~,~,~, ~, ~, ~, f,thetaX] = ...
                            PLI_pair([iData, dataToUse_filt(:,jSensor)], Fs, freqBands{iF},thetaX);                                                                

                    end

                end
                % run psd
                tmp= GC_psd(iData, nfft, Fs, hamming(windSampleSize));
                if isempty(ppsd)                
                    ppsd(nSensors,length(tmp),nWindows, nFreqs)=0;
                end;
                ppsd(iSensor, :, iW,iF) = tmp;

            end

            % symmetrise matrices
            %PLVw(:,:,iW,iF)          = PLVw(:,:,iW)          + PLVw(:,:,iW).';


            clear dataToUse_filt;
        end; % nFreqs
    end % nWindows

    PSD=zeros(nFreqs,nSensors);
    PLV=zeros(nFreqs,nSensors,nSensors);

    for iFtop=1:nFreqs,
        freqBand=freqBands{iFtop};
        % frequencies of relevance

        f(:,2:end) = []; % remove duplicate data required by parfor
        fIn = f<=freqBand(2) & f >= freqBand(1);

        %% PSD
        PSD(iFtop,:)    = mean(mean(ppsd(:,fIn,:,iFtop),3),2);

        %% Phase measures
        PLV(iFtop,:,:)    = abs(mean(PLVw(:,:,:,iFtop),3));

    end %nFreqs

else
    PSD=nan(nFreqs,nSensors);
    PLV=nan(nFreqs,nSensors,nSensors);
end

end%PLI_band2

function [PLVw, sindTheta, signsindTheta, abssindTheta, cpsd, MI, PLT, f, thetaX] = PLI_pair(data, Fs, freqBand, thetaX)
[nSamples, checkMe] = size(data);
assert(2 == checkMe);
nfft = nSamples;
cpsd=[];
signsindTheta=[];
sindTheta=[];
abssindTheta=[];

% Hilbert transform
if isempty(thetaX),
    htx = ROInets.demean(hilbert(data(:,1), nfft));
    thetaX = angle(htx(1:ROInets.rows(data)));
    thetaX(1:10) = []; thetaX(end-9:end) = [];

end
hty = ROInets.demean(hilbert(data(:,2), nfft));
thetaY = angle(hty(1:ROInets.rows(data)));
% remove edge effects
thetaY(1:10) = []; thetaY(end-9:end) = [];
% phase
relativePhase = thetaX - thetaY;
%phiX = mod(thetaX, 2 * pi) - pi;
%phiY = mod(thetaY, 2 * pi) - pi;
% entropy of phase
%I    = mutualinfo(phiX, phiY);
%H    = jointentropy(phiX,phiY);
I=[];
H=[];

% psd
%window = hamming(nSamples);
f     = frequencygrid('half', nfft, Fs, false);


% phase difference between -pi and pi
%dTheta = mod(relativePhase, 2 * pi) - pi;
%sindTheta = mean(sin(dTheta));
%signsindTheta = mean(sign(sin(dTheta)));
%abssindTheta  = mean(abs(sin(dTheta)));
PLVw          = mean(exp(1i .* relativePhase(:)));

% average over frequencies
%fIn  = f<=freqBand(2) & f >= freqBand(1);
% cpsd = cpsd(fIn);
% f    = f(fIn);

% mutual information of phases, as outlined in eq 10 of Palus 1997
% "Detetecting phase synchronization in noisy systems" 
% and using normalisation from Wilmer 2012 "Time-Delayed Mutual information
% of the phase as a measure of functional connectivity"
%MI = (I./H);
MI=[];

if(0)
    % from here on, compute PLT
    % step 4, is this greater or smaller than zero?
    dTheta2=sign(sin(dTheta));
    % step 5, find crossings
    cross = zeros(numel(dTheta2),1);
    for b=1:numel(dTheta2)-1,
        if dTheta2(b) ~=dTheta2(b+1),
            cross(b+1) = 1;
        end%if
    end%if
    % step 6, compute time between crossings
    crossings = find(cross);
    time      = diff(crossings)./Fs;
    % step 7, if time between crossings too small, no synchronization
    meanFreq    = mean(freqBand);
    thresh      = repmat(1.0/meanFreq, numel(time), 1);
    if ~isempty(time),
    index       = sign(time-thresh) < 1;
    time(index) = [];
    end
    if ~isempty(time),
        PLT = mean(1-exp(-time));
    else
        PLT = 0;
    end%if
else
    PLT=[];
end;

end%PLI_pair

function windInd = get_windows(nSamples, windowSize, Fs)
%get_windows returns logical indices pulling out overlapping windows from a datset

windowLength = ceil(windowSize * Fs);

% desire 0% overlap.
nWindows = floor( nSamples / windowLength );

windInd = false(nSamples, nWindows);
for iW = 1:nWindows,
    windInd((iW-1) * windowLength  + 1 : min((iW-1) * windowLength + windowLength, nSamples),iW) = true;
end
end

function N = vector_norm(M)
%vector_norm takes the norms of columns of matrix method

for iCol = cols(M):-1:1,
    N(iCol) = sqrt(M(:,iCol)' * M(:,iCol));
end
end





function [ PXY ] = GC_cpsd(x,y,nfft,Fs, window)
% CPSD for pre-chunked data.

% window the data
xw = bsxfun(@times, x, window);
yw = bsxfun(@times, y, window);

U = window'*window; % compensates for power of window

% Compute the periodogram power spectrum [Power] estimate
% A 1/N factor has been omitted since it cancels
assert(nfft>=ROInets.rows(x));

Xx = fft(xw,nfft);
Yy = fft(yw,nfft); 

% P = Xx.*conj(Xx)/U;      % Auto spectrum.
SXY = bsxfun(@times,Xx,conj(Yy))/U;  % Cross spectrum mean square power.

% make one-sided and convert to psd
% rip from computepsd

% Generate the one-sided spectrum [Power] if so wanted
if rem(nfft,2),
    select = 1:(nfft+1)/2;  % ODD
    Sxy_unscaled = SXY(select,:); % Take only [0,pi] or [0,pi)
    Sxy = [Sxy_unscaled(1,:); 2*Sxy_unscaled(2:end,:)];  % Only DC is a unique point and doesn't get doubled
else
    select = 1:nfft/2+1;    % EVEN
    Sxy_unscaled = SXY(select,:); % Take only [0,pi] or [0,pi)
    Sxy = [Sxy_unscaled(1,:); 2*Sxy_unscaled(2:end-1,:); Sxy_unscaled(end,:)]; % Don't double unique Nyquist point
end

% Compute the PSD [Power/freq]
PXY = Sxy./Fs; % Scale by the sampling frequency to obtain the psd
end

function [ PXX ] = GC_psd(x, nfft, Fs, window)
% PSD for pre-chunked data.

% window the data
xw = bsxfun(@times, x, window);

U = window'*window; % compensates for power of window

% Compute the periodogram power spectrum [Power] estimate
% A 1/N factor has been omitted since it cancels
assert(nfft>=ROInets.rows(x));

Xx = fft(xw,nfft);

SXX = Xx.*conj(Xx)/U;      % Auto spectrum.

% make one-sided and convert to psd
% rip from computepsd

% Generate the one-sided spectrum [Power] if so wanted
if rem(nfft,2),
    select = 1:(nfft+1)/2;  % ODD
    Sxx_unscaled = SXX(select,:); % Take only [0,pi] or [0,pi)
    Sxx = [Sxx_unscaled(1,:); 2*Sxx_unscaled(2:end,:)];  % Only DC is a unique point and doesn't get doubled
else
    select = 1:nfft/2+1;    % EVEN
    Sxx_unscaled = SXX(select,:); % Take only [0,pi] or [0,pi)
    Sxx = [Sxx_unscaled(1,:); 2*Sxx_unscaled(2:end-1,:); Sxx_unscaled(end,:)]; % Don't double unique Nyquist point
end

% Compute the PSD [Power/freq]
PXX = Sxx./Fs; % Scale by the sampling frequency to obtain the psd
end







%--------------------------------------------------------------------------
function w = frequencygrid(Range,Npts,Fs,CenterDC)
% Compute the frequency grid.

% Compute the whole frequency range, e.g., [0,2pi) to avoid round off errors.
if isempty(Fs),
    Fs = 2*pi;
end
freq_res = Fs/Npts;
w = freq_res*(0:Npts-1);

% There can still be some minor round off errors in the frequency grid.  
% Fix the known points, i.e., those near pi and 2pi.
Nyq = Fs/2;
half_res = freq_res/2; % half the resolution

% Determine if Npts is odd and calculate half and quarter of Npts.
[isNPTSodd,halfNPTS,ishalfNPTSodd,quarterNPTS] = NPTSinfo(Npts);

if isNPTSodd,
    % Adjust points on either side of Nyquist.
    w(halfNPTS)   = Nyq - half_res;
    w(halfNPTS+1) = Nyq + half_res;
else
    % Make sure we hit Nyquist exactly, i.e., pi or Fs/2 
    w(halfNPTS) = Nyq;
end
w(Npts) = Fs-freq_res;

% Get the right grid based on range, centerdc, etc.
w = finalgrid(w,Npts,Nyq,Range,CenterDC,isNPTSodd,ishalfNPTSodd,halfNPTS,quarterNPTS);
end
%--------------------------------------------------------------------------
function [isNPTSodd,halfNPTS,ishalfNPTSodd,quarterNPTS] = NPTSinfo(NPTS)
% Determine if we're dealing with even or odd lengths of NPTS, 1/2 NPTS,
% and 1/4 NPTS.

% Determine if Npts is odd.
isNPTSodd = false;
if rem(NPTS,2),
    isNPTSodd = true;
end

% Determine half the number of points.
if isNPTSodd,   halfNPTS = (NPTS+1)/2;  % ODD
else            halfNPTS = (NPTS/2)+1;  % EVEN
end

% Determine if half Npts is odd.
ishalfNPTSodd = false;     
if rem(halfNPTS,2),        
    ishalfNPTSodd = true;  
end

% Determine a quarter of the number of points.
if ishalfNPTSodd,  quarterNPTS = (halfNPTS+1)/2;  % ODD
else               quarterNPTS = (halfNPTS/2)+1;  % EVEN
end
end
%--------------------------------------------------------------------------
function w = finalgrid(w,Npts,Nyq,Range,CenterDC,isNPTSodd,ishalfNPTSodd,halfNPTS,quarterNPTS)
% Calculate the correct grid based on user specified values for range,
% centerdc, etc.

switch lower(Range)
    case 'whole',
        % Calculated by default.% [0, 2pi)

        if CenterDC,          % (-pi, pi] even or (-pi, pi) odd
            if isNPTSodd,  negEndPt = halfNPTS;
            else           negEndPt = halfNPTS-1;
            end
            w = [-fliplr(w(2:negEndPt)), w(1:halfNPTS)];
        end
        
    case 'half'            
        w = w(1:halfNPTS);      % [0, pi] even or [0, pi) odd
        
        % For even number of points that are not divisible by 4 you get
        % less one point to avoid going outside the [-pi/2 pi/2] range.
        if CenterDC,            % [-pi/2, pi/2] even (-pi/2, pi/2) odd 
            if ishalfNPTSodd,
                negEndPt = quarterNPTS;
            else
                quarterNPTS = quarterNPTS-1; % Avoid going over pi/2
                negEndPt = quarterNPTS;
            end
            w = [-fliplr(w(2:negEndPt)), w(1:quarterNPTS)];
            if ~rem(Npts,4),
                % Make sure we hit pi/2 exactly when Npts is divisible
                % by 4! In this case it's due to roundoff.
                w(end) = Nyq/2;
            end
        end
    otherwise
        error(message('signal:psdfreqvec:InternalError'));
end
w = w(:);  % Return a column vector.
end

% [EOF]
