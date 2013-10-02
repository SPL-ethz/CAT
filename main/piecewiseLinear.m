function [y] = piecewiseLinear(varargin)
% % piecewiseLinear(tpoints,ypoints,t) or piecewiseLinear([tpoints;ypoints],y)
% evaluates a piecewise linear function with nodes at (tpoints,ypoints) at
% time(s) t


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
        'Your PWC profile has an unsorted time');
end
    
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

y = a(I)+m(I).*t;
y = reshape(y,size(t));