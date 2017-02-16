function [rsq, p] = rsq(x,y)
% r^2 between x and y (ignores NaN)

% [tmp, p, ~, ~, Tstat] = corrcoef(x,y,'rows','complete');
[tmp, p] = corrcoef(x,y,'rows','complete');
rsq = tmp(2).^2;
p = p(2);
