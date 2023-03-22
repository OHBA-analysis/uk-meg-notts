# uk-meg-notts

Scripts to fit an HMM to the uk_meg_notts dataset using the MATLAB OSL (with HMM-MAR) and the python OSL (with osl-dynamics).

Repos:

- MATLAB OSL: https://ohba-analysis.github.io/osl-docs
- HMM-MAR: https://github.com/OHBA-analysis/HMM-MAR
- Python OSL: https://github.com/OHBA-analysis/osl
- osl-dynamics: https://github.com/OHBA-analysis/osl-dynamics

Full matlab pipeline:
```
matlab -nodesktop
>> osl_startup;
>> convert_raw_data;
>> preprocess_data;
>> source_reconstruct;
>> prepare_te_pca_data;
>> fit_hmm;
>> analyse_hmm_fit;
>> save_maps;
```
Note, the MATLAB OSL package is obselete and is not being maintained.

Fully python pipeline (recommended):
```
>> python preprocess.py
>> python source_reconstuct.py
>> python sign_flip.py
>> python train_hmm.py
>> python calc_multitaper.py
>> python plot_networks.py
```
