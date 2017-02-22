function data = mapVal(data, varargin)
% data = map(data, dict)
% dict is 2 x N; (source-->target)
%  OR
% data = map(data, source, target)

oldData = data;
if nargin==2
   source = varargin{1}(1,:);
   target = varargin{1}(2,:);
end

if nargin==3
   source = varargin{1};
   target = varargin{2};
end

for dat = 1:size(source,2)
   data(oldData==source(dat)) = target(dat);
end
