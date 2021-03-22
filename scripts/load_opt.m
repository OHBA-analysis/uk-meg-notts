%
% Loads data from the OSL preprocessing tool
%
load([dirs.vars '/optName'], 'optName');

opt = osl_load_opt([dirs.data, '/' optName]);
opt = opt_consolidate_results(opt);

clear optName;
