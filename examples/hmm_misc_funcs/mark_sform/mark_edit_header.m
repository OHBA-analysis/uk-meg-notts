function [in,out] = mark_edit_header( input_file, output_file, pixdim, offset, qfac )

    if nargin < 5, qfac=-1; end
    if nargin < 4, offset=zeros(1,3); end
    if nargin < 3, pixdim=ones(1,3); end

    assert( ischar(input_file) && ischar(output_file), 'First and second inputs must filenames.' );
    assert( isnumeric(pixdim) && numel(pixdim)==3 && all(pixdim>0), 'Third input should be a 1x3 positive vector of pixdims.' );
    assert( isnumeric(offset) && numel(offset)==3, 'Fourth input should be a 1x3 translation vector.' );
    assert( isnumeric(qfac) && isscalar(qfac), 'Fifth input should be a scalar q-factor.' );

    % load input nifti file
    in = load_untouch_nii(input_file);
    
    % information needed for solving the orientation
    dat = in.hdr.hist;
    dim = in.hdr.dime;
    
    % for affine transformations only, the simplest is to use the s-form
    dat.qform_code = 0;
    dat.sform_code = 4;
    %dat.magic = 'n+1';
    
    % set pixdims
    dim.pixdim(1) = qfac;
    dim.pixdim([2 3 4]) = pixdim;
    
    % set s-form
    dat.srow_x = [ qfac*pixdim(1)         0         0 offset(1) ];
    dat.srow_y = [              0 pixdim(2)         0 offset(2) ];
    dat.srow_z = [              0         0 pixdim(3) offset(3) ];
    
    % set corresponding q-form too just to be sure
    dat.qoffset_x = offset(1);
    dat.qoffset_y = offset(2);
    dat.qoffset_z = offset(3);
    
    dat.quatern_b = 0;
    dat.quatern_c = 0;
    dat.quatern_d = 0;
    
    % save output file with modified headers
    out = in;
    out.hdr.hist = dat;
    out.hdr.dime = dim;
    
    fprintf('Saving output file: "%s"\n',output_file);
    save_untouch_nii( out, output_file );
    
end