function res = allvaleq ( A , tol )

%ALLVALEQ   Check if all values in a matrix are equal
%   ALLVALEQ( A ) returns logical true if all values in A are equal (within
%   a tolerance) and logical false otherwise. By default, the tolerance is
%   1e-10 . The tolerance is multiplied by max(abs(A)).
%
%   ALLVALEQ( A , tol ) uses the tolerance tol instead of the default value

if nargin < 2
    tol = 1e-10;
end % if

res = all(abs( A(:) - mean( A(:) ) ) <= tol*max(abs(A(:))) );

end % function