%
% Fits an HMM
%
% Must run run_preprocessing.m, run_source_reconstruction.m, run_data_preparation.m
% before this script.
%
setup;
get_initial_hmm_covariances;
fit_hmm;
compute_spectra;
spectral_factorisation;
get_spatial_maps;
