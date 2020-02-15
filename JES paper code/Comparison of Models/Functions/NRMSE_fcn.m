%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to compute the NRSME between two variables
%
% On Trade-offs between Computational Complexity and Accuracy of
% Electrochemistry-based Battery models
%
% Authors: Z. Khalik, H.J. Bergveld, M.C.F. Donkers
%
% This file is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function NRMSE = NRMSE_fcn(x,y)
N = length(x); 
NRMSE = sqrt(1/N*sum((x-y).^2))/(max(0.5*(x+y))-min(0.5*(x+y))); 
end