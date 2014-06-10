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

if length(varargin) >= nargin(func)
    
    % Evaluate the function with the first n arguments
    out = feval(func,varargin{1:nargin(func)});
    
else
    error('CAT:evalanonfunc:notenoughinput',...
        'Not enough input arguments')
end % if else

end % function