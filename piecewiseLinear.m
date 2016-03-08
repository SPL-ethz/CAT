function [y] = piecewiseLinear(varargin)
% % piecewiseLinear(tpoints,ypoints,t) or piecewiseLinear([tpoints;ypoints],y)
% evaluates a piecewise linear function with nodes at (tpoints,ypoints) at
% time(s) t

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



if nargin == 3
    tpoints = varargin{1}; tpoints = tpoints(:);
    ypoints = varargin{2}; ypoints = ypoints(:);
    t = varargin{3};
elseif nargin == 2
    tpoints = varargin{1}(1,:); tpoints = tpoints(:);
    ypoints = varargin{1}(2,:); ypoints = ypoints(:);
    t = varargin{2};
end

if ~issorted(tpoints)
    warning('piecewiseLinear:unsortedTime',...
        'Your PWL profile has an unsorted time');
end

t(t<tpoints(1) & t>=tpoints(1)-1e-5) = tpoints(1);

    
tpoints = tpoints(:);
ypoints = ypoints(:);
t = t(:);

% find slopes between points
m = diff(ypoints)./diff(tpoints);

% find abscissa
a = ypoints(1:end-1)-m.*tpoints(1:end-1);

for i = 1:length(t)
    I(i)=find(t(i)>=tpoints,1,'last');
end

I(I==length(ypoints)) = length(ypoints)-1;
a = a(I);m = m(I);
y = a(:)+m(:).*t;
y = reshape(y,size(t));