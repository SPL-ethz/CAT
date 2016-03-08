function out = evalanonfunc(func,varargin)

% evalanonfunc
%
% Evaluate an anynomous function 'func' using the input arguments given.
%
% evalanonfunc checks the number of inputs which the specified anonymous
% function accepts, and passes the first n of those inputs to the function,
% returning the output.
%
% SEE ALSO
% eval

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


if length(varargin) >= nargin(func)
    
    % Evaluate the function with the first n arguments
    out = feval(func,varargin{1:nargin(func)});
    
else
    error('CAT:evalanonfunc:notenoughinput',...
        'Not enough input arguments')
end % if else

end % function