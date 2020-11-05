function [] = add_specs(in_fname, out_fname, force, ...
                        mt_opt, ga_opt, mar_opt, gg_opt)
%ADD_SPECS: add spectral information to HMM-MAR structure
addpath(genpath(strcat(getenv('HOME'), '/HM/SW/HMM-MAR')));

% Check that HMMMMAR output has all the required metadata as per latest
% standard, and add spectral magnitudes if needed, according to parameters.
%
% Adds: multitaper spectra, MAR spectra and granger causality
%
% It makes sense to call this function in a loop when in is just one
% file (no trials), when there are several in files, we must calculate
% how they break up into trials
%
%
% INPUT
%
% in_fname      filename of .h5 file containing the data
% out_fname     output filename (.mat)
%
% force         overwrite all existing spectral MT, MAR, Granger info
%
% mt_opt        multitaper spectral options
%     .Fs            sampling rate of data
%     .tapers        [TW prod, #tapers], generally #tapers=2*TW-1
%     .win           windows for the multitaper spectra
%     .fpass         spectral domain (e.g. [0 280]), in Hz
%     .p             confidence level for the uncertainty in spectra
%     .numIterations for the Wilson iterative factorization
%     .tol           for the Wilson iterative factorization
%
% ga_opt multitaper opts for the responsibilities (gamma, ga). See
%        mt_opt. The spectra of the gammas will be computed as if they
%        were separate channels, with no weigths
%        (ga_opt.Gamma=ones(...)). Option ga_opt.to_do is a 2-vector
%        indicating whether or not (1, 0 resp) to compute coherences
%        and PDCs
%
%
% mar_opt        mar (parametric) spectral options
%     .Nf           number of frequencies in fpass to sample for the MAR PSD
%     .fpass        spectral domain for MAR spectra
% gg_opt         Granger options
%     .gg_alpha     alpha confidence level for Granger summary measure
%
% OUTPUT
% metadata mt_opt, ga_opt, mar_opt, gg_opt and output
% spec_mt, spec_ga, spec_mar and granger appended at out_filename

% Álvaro Tejero-Cantero, tejero at bio.lmu.de (2015)


%    XXX T could be passed from main, but for the moment we get it here
% so as to change interfaces the least
% when trial information is available, T should be overwritten
% XXX is this really needed, given that T is loaded from the out-fname
% below?
%  X = h5read(report.in_fname, '/X');
%
% try
%         T = h5read(report.in_fname, '/T');
% catch err
%         [T, ndim] = size(X, 1);
%         fprintf('No trial information found in %s/T', report.in_fname);
%         fprintf('\n');
%         fprintf('Will assume recording continuous of length %d', T);
% end;

load(out_fname); % need T, Gamma (for mt) and hmm (for mar) from here

X = h5read(in_fname, '/X');


%%% add multitaper spectra
if force || ~exist('spec_mt', 'var')
    disp('No MT signal spectra found or force=1.');

    mt_opt.Gamma = Gamma; % weight by occupancy
    spec_mt = hmmspectramt(X, T, mt_opt);

    mt_opt = rmfield(mt_opt, 'Gamma'); %spare saving it
    save(out_fname, 'mt_opt', 'spec_mt', '-append');
end


%%% add spectra of responsibilities

if force || ~exist('spec_ga', 'var')
    disp('No MT Gamma spectra found or force=1.');

    Tcut = T;
    % figure out furthest regressor from current sample
    real_order = formorders(hmm.train.order, hmm.train.orderoffset, ...
                            hmm.train.timelag, hmm.train.exptimelag);

    for in=1:length(T); Tcut(in) = Tcut(in) - real_order(end); end;

    % spectra of Gamma with no weights: indeed, ga_opt.Gamma=ones(...) by
    % default; hijacks normal operation where first argument is signal
    % and ga_opt.Gamma contains the weights.
    spec_ga = hmmspectramt(Gamma, Tcut, ga_opt);

    save(out_fname, 'ga_opt', 'spec_ga', '-append');
end


%%%% add MAR spectra
if force || ~exist('spec_mar', 'var')
    disp('No MAR spectra found or force=1.');

    mar_opt.Gamma = Gamma; % needed since Jan 2015
    spec_mar = hmmspectramar(hmm, X, T, mar_opt);

    mar_opt = rmfield(mar_opt, 'Gamma'); %spare saving it
    save(out_fname, 'mar_opt', 'spec_mar', '-append');
end

end % function



%%%% add Granger causality - INACTIVE - GRANGER IN ALPHA
% if force || ~exist('granger','var')
%     disp('No Granger metrics found of force=1.');
%     granger = struct();
%     [granger.pval,granger.sig,granger.F,...
%         granger.f] = tgrangerc (X,T,Gamma, 100,gg_opt.ggalpha,...
%                                 hmm.train.Fs,1,0)
%     granger.ggalpha = gg_opt.alpha;
%     save(out_fname, 'granger', '-append');
% end
