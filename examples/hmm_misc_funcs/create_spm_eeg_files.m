function [ file_names_out ] = create_spm_eeg_files( S )

% [ file_names_out ] = create_spm_eeg_files( S )
%
% S.datasets % cell array of 2D matrices of data (space x samples)
% S.res: temp resolution of data
% S.template_spm_filename
% S.new_filenames % cell array of output file names

%%%%%%%%%%
%% extract each subjects data and put it into a clone of a dummy spm file
file_names_out=S.new_filenames;

dummy_spm_file=spm_eeg_load(S.template_spm_filename);

for ss=1:length(S.datasets)    
    Sc=[];
    Sc.D=dummy_spm_file;
    Sc.newdata=S.datasets{ss}(:,:,:);
    Sc.newname=[S.new_filenames{ss} '.dat'];
    Sc.time=S.res:S.res:S.res*size(S.datasets{ss},2);
    Sc.chantype='MEG';
    osl_change_spm_eeg_data( Sc );
        
end

end

