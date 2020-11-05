%
% Full Pipeline
%
setup;

%
% Preprocessing
%
get_subjects;
convert_ds_files;
convert_mri_files;
run_opt;
split_into_sessions;

%
% Source reconstruction
%
run_beamforming;
run_parcellation;
fix_dipole_sign;
save_preprocessed_data;

%
% Prepare data for the HMM (time embedding and PCA)
%
prepare_data;

%
% HMM analysis
%
get_initial_hmm_covariances;
fit_hmm;
compute_spectra;
spectral_factorisation;
get_spatial_maps;
compute_state_statistics;
