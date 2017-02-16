% clean up workspace: command window, figures, variables
try
   if ~isempty(findall(0,'Type','Figure'))
      clf
   end
end
clear
clc
clear global