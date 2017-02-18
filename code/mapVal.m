function data = mapVal(data, varargin)
% data = map(data, dict)
% dict is 2 x N; (source-->target)
%  OR
% data = map(data, source, target)

%function data = map(data, varargin)

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



% oldData = data;
% uniData = unique(data);
%
% for dat = 1:length(uniData)
%    idx = oldData==uniData(dat);
%    dictIdx = find(dict(1,:)==uniData(dat));
%    if ~isempty(dictIdx)
%       data(idx) = dict(2,idx);
%    else
%       data(idx) = nan;
%    end
% end
