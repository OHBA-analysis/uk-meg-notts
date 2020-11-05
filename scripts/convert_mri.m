function [] = convert_mri(mriFileName, niftiFileName)
%CONVERT_MRI  Converts .mri files to nifti
%
% convert_mri(mriFileName, niftiFileName) converts CTF .mri file to nifti
%   format. 
%
%   Conversion courtesy George O'Neill at Notthingham. 

%	Copyright 2014 OHBA
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


%	$LastChangedBy: GilesColclough $
%	$Revision: 763 $
%	$LastChangedDate: 2015-10-21 11:52:19 +0100 (Wed, 21 Oct 2015) $
%	Contact: giles.colclough@magd.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 05-Nov-2014 12:35:40

% input checking
narginchk(2,2);
assert(ischar(mriFileName) && exist(mriFileName, 'file'), ...
       [mfilename ':InputFileNotExist'], ...
       'Can''t find input file mriFileName. \n');
assert(ischar(niftiFileName), ...
       [mfilename ':BadNiftiString'], ...
       'niftiFileName must be a character string. \n');
[fileDir, fileName, fileExt] = fileparts(niftiFileName);
assert(isempty(fileDir) || isdir(fileDir), ...
       [mfilename ':TargetDirNotExist'], ...
       'Target directory for nifti file does not exist. \n');
assert(strcmpi(fileExt, '.nii') || (strcmpi(fileExt, '.gz') && strcmpi(fileName(end-3:end), '.nii')), ...
       [mfilename ':BadNiftiExt'], ...
       'niftiFileName must have a nifti extension. \n');

% Begin conversion
mri = ft_read_mri(mriFileName); 

mri.anatomy = flip(mri.anatomy,1);

mri.anatomy = flip(mri.anatomy,2);

mri.anatomy = flip(mri.anatomy,3);

mri.anatomy = flip(mri.anatomy,1);

 
% write out
ft_write_mri(niftiFileName, mri.anatomy, 'dataformat','nifti');

% This bit added by MWW to make sure hdr info is correct
nii = load_untouch_nii(niftiFileName);

[pth nm ext]=fileparts(niftiFileName);

% We assume qform and sform are the same, and correspond to the
% transform from voxel coords to native scanner space.
% For some reason rhino does not cope well if native scanner
% space coords are too different from MNI. So to help with this
% we change the offset to get the native scanner space coords
% looking a bit more like MNI (note that the chosen native 
% scanner space coordinate system is arbitrary wrt what rhino 
% needs to do). (presumably it is a local minima
% issue in rhino - and ideally should be fixed in there).
% We need to make sure the qformcode and sformcode is 1.
% As the output from edit_header
% sets an invalid code which rhino will not like.
% See convert_mri for more on this
offset=[128 -128 -90];
[in,out] = edit_header( niftiFileName, fullfile(pth,nm), nii.hdr.dime.pixdim(2:4), offset, -1 );

%%%%%%% 
% added by MWW Aug 2019
%%%%%%% 
% add correct sform code.
% the output from edit_header has an invalid qform and
% seemingly correct sform (we can check sform is working as
% when it is viewed in fsleyes, all of the labels seem to be
% correct (although can not be sure about left-right or the offset) ).
% However, the sform code is set to 4. For rhino to be happy
% the sform code needs to be set to 1. So we manually enforce
% this now. (Note that this call also seems to set the qform to be
% the same as the sform):
% fslorient -setsformcode 1 imagename
runcmd(['fslorient -setsformcode 1 ' niftiFileName]);

%runcmd(['rm -f ' fullfile(pth,[nm,'.nii'])]);
%runcmd(['gunzip ' fullfile(pth,[nm,'.nii.gz'])]);
            

end%convert_mri
% [EOF]
