function [rsq, p] = rsq(x,y)
% [rsq, p] = rsq(x,y)
% r^2 between x and y (ignores NaN)

[tmp, p] = corrcoef(x,y,'rows','complete');
rsq = tmp(2).^2;
p = p(2);
