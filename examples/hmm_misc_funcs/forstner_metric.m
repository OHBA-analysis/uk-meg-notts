function d = forstner_metric(A, B)
%FORSTNER_METRIC
% Forstner's metric for the difference between two covariance matrices
% 
%   D = FORSTNER_METRIC(A, B) calculates Forstner's metric for the 
%   difference between two positive definite covariance matrices. 
%    
%   For the derivation, see Förstner, W. & Moonen, B. A metric for 
%   covariance matrices. Quo vadis geodesia 113–128 (1999).
%   http://www.uni-stuttgart.de/gi/research/schriftenreihe/quo_vadis/pdf/foerstner.pdf

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
%	$Revision: 1 $
%	$LastChangedDate: 2013-07-22 14:22:00 +0100 (Mon, 22 Jul 2013) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 22-Jul-2013 14:22:00

%% Check inputs
isMatrixInput = @(X) isnumeric(X) && ismatrix(X);
if ~isMatrixInput(A) || ~isMatrixInput(B),
    error([mfilename ':NotMatrixInput'], ...
          'Arguments must be matrices\n');
end%if

if ~isequal(size(A), size(B)),
    error([mfilename ':MatricesDifferInSize'], ...
          'Input matrices must be of equal size\n');
end%if

[n, m] = size(B);
if ~isequal(n, m),
    error([mfilename ':NotSquareInput'], ...
          'Arguments must be square matrices\n');
end%if

%if ~isposdef(A) || ~isposdef(B)
%    error([mfilename ':NotPosDef'], ...
%          'Arguments must be positive definite matrices\n');
%end%if

%% Calculate metric
E = eig(B, A, 'chol');

d = sqrt(sum(log(E).^2));

end%forstner_metric
% [EOF]