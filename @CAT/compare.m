function [results,same,different] = compare(CAT1,CAT2)
%% [results,same,different] = CAT.compare(CAT1,CAT2)
% compare is a (static) method of the CAT class. It compares two CAT
% objects and indicates similarities and differences between the two
% objects.
% If no output is assigned, compare prints an analysis of the comparison in
% the command window. If outputs are assigned, compare returns a logical
% vector 'results' whose elements are true for identical properties and
% false otherwise (order is the same as in properties(CAT)), a cell array
% 'same' containing the names of properties that are identical, a cell
% array 'different' containing the names of properties that are different.

% Copyright 2012-2014 David Ochsenbein, Martin Iggland; 2014-2016 David Ochsenbein
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


p = properties(CAT1);
fprintf('Property-->Result \n-------------------\n')
result = zeros(length(p),1);
for i = 1:length(p)
    
    if isa(CAT1.(p{i}),'function_handle') && isa(CAT2.(p{i}),'function_handle')
        str1 = func2str(CAT1.(p{i}));
        str2 = func2str(CAT1.(p{i}));
    elseif isa(CAT1.(p{i}),'Distribution') && isa(CAT2.(p{i}),'Distribution')
        str1 = data2str(CAT1.(p{i}));
        str2 = data2str(CAT1.(p{i}));
    end
        result(i) = isequal(str1,str2);
    
    
    if nargout == 0
        if result(i)
            resultstr = 'EQUAL';
        else
            resultstr = 'DIFFERENT';
        end

        fprintf(strcat(p{i},'-->',resultstr,'\n'))
    end
    
end

if nargout>0
    results = result;
    issame = 1;isdiff = 1;
    same = [];different = [];
    for i = 1:length(result)
        if result(i)
            same{issame} = p{i};
            issame = issame +1;
        else
            different{isdiff} = p{i};
            isdiff = isdiff + 1;
        end
    end
end
