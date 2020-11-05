function VV = GC_normalise_vectors(V, dim)
%GC_NORMALISE_VECTORS normalises rows or columns of a matrix
%
% NORMV = GC_NORMALISE_VECTORS(V, DIM) normalises vectors along dimension 
%   DIM of V. If DIM=1, this treats columns as vectors and if DIM=2, this
%   treats rows as vectors. 
%
% NORMV = GC_NORMALISE_VECTORS(V) is the same as GC_NORMALISE_VECTORS(V, 1)
%
% Example:
%   If X = [3 3 3]
%   Then GC_NORMALISE_VECTORS(X,1) is [1 1 1] and GC_NORMALISE_VECTORS(X,2)
%   is [0.5774 0.5774 0.5774]. 


%	Copyright 2013 Giles Colclough
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
%	$Revision: $
%	$LastChangedDate: $
%	Contact: giles.colclough 'at' eng.ox.ac.uk
%	Originally written on: GLNXA64 by Giles Colclough, 25-Nov-2013 11:49:41

if nargin < 2 || ~exist('dim', 'var') || isempty(dim),
    dim = 1;
end%if

VV = bsxfun(@rdivide, V, sqrt(sum(V.^2, dim) ./ size(V, dim)));
end%GC_normalise_vectors