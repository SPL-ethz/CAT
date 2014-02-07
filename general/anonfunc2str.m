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

% Get string of function
str = func2str(func);

% Use regexp to get rid of @(...) part at beginning
S = regexp(str,'^(@\(.*\))(.*)','tokens');

if ~isempty(S{1})

str = S{1}{2};
argsin = S{1}{1};

else
    error('anonfunc2str:extract:nomatch',...
        'The function could not be translated. Is it an anonymous function?');
end % if else

end % function anonfunc2str