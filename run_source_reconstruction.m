%
% Source Reconstruction
%
% Must run run_preprocessing.m before this script.
%
setup;
run_beamforming;
run_parcellation;
fix_dipole_sign;
save_preprocessed_data;
