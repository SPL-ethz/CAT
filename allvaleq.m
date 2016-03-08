function res = allvaleq ( A , tol )

%ALLVALEQ   Check if all values in a matrix are equal
%   ALLVALEQ( A ) returns logical true if all values in A are equal (within
%   a tolerance) and logical false otherwise. By default, the tolerance is
%   1e-10 . The tolerance is multiplied by max(abs(A)).
%
%   ALLVALEQ( A , tol ) uses the tolerance tol instead of the default value

% Copyright 2015-2016 David Ochsenbein
% Copyright 2012-2014 David Ochsenbein, Martin Iggland
% 
% This file is part of CAT.
% 
% CAT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
% 
% CAT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


if nargin < 2
    tol = 1e-10;
end % if

res = all(abs( A(:) - mean( A(:) ) ) <= tol*max(abs(A(:))) );

end % function