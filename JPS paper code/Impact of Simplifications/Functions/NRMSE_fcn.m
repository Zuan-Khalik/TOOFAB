%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Model Simplifications and Its Impact on Computational Complexity for an 
% Electrochemistry-Based Battery Modeling Toolbox
%
% Authors: Z. Khalik, M.C.F. Donkers, H.J. Bergveld
%
% This file is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function NRMSE = NRMSE_fcn(x,y)
N = length(x); 
NRMSE = sqrt(1/N*sum((x-y).^2))/(max(0.5*(x+y))-min(0.5*(x+y))); 
end