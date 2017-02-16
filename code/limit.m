function x = limit(x, varargin)
%x = limit(x, [low=0], [up=1])
if nargin==1
   up = 1;
   low = 0;
else
   low = varargin{1};
   up = varargin{2};
end
x = max(x,low); % Bound elements from below, x >= lowerBound
x = min(x,up); % Bound elements from above, x <= upperBound 
