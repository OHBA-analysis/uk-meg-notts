# uk-meg-notts

Scripts to fit an HMM to the uk_meg_notts dataset using OSL.

Full pipeline:
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
