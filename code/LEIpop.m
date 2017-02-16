function [er, param] = LEIpop(param, pa)
% [er, param] = LEIpop(param, pa)


% rescale parameters from [-1 1] to [lb ub]
for i = 1:size(param,2)
   param(:,i) = (param(:,i)+1)/2*(pa.ub(i) - pa.lb(i)) + pa.lb(i);
end

% get MSE for each individual from objFunind
inds = size(param,1);
er = zeros(1,inds);
for ind = 1:inds
   er(ind) = pa.objFunInd(param(ind,:), pa);
end
% convert error to a normalized fitness measure
er = 1-er/var(pa.meanResp);
