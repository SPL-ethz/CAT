function [str,argsin] = anonfunc2str(func)

% anonfunc2str
%
% str = anonfunc2str(anonfunction);
% returns a string representation of an anonymous function, without the
% input argument definition! (Difference to func2str)
%
% [str,argsin] = anonfunc2str(anonfunction);
% returns the string for the input arguments definition as well
%
% SEE ALSO
% func2str, data2str

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


% Get string of function
str = func2str(func);

% Use regexp to get rid of @(...) part at beginning
S = regexp(str,'^(@\([\w,]*?\))(.*)','tokens');

if ~isempty(S{1})
    
    str = S{1}{2};
    argsin = S{1}{1};
    
else
    error('anonfunc2str:extract:nomatch',...
        'The function could not be translated. Is it an anonymous function?');
end % if else

end % function anonfunc2str