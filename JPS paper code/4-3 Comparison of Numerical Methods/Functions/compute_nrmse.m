%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to compute the NRMSE values of the output and several internal
% states between two simulation output files
%
% Model Simplifications and Its Impact on Computational Complexity for an 
% Electrochemistry-Based Battery Modeling Toolbox
%
% Authors: Z. Khalik, M.C.F. Donkers, H.J. Bergveld
%
% This file is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function NRMSE = compute_nrmse(out1,out2)
nt = length(out1.V); 
nx = out1.param.grid.nn+out1.param.grid.ns+out1.param.grid.np;
nn = out1.param.grid.nn; 
nnp = out1.param.grid.nn+out1.param.grid.np;

NRMSE.V = NRMSE_fcn(out1.V,out2.V)*1000; 

phie1 = reshape(out1.phie, [1 nx*nt]); 
phie2 = reshape(out2.phie, [1 nx*nt]); 
NRMSE.phie = NRMSE_fcn(phie1,phie2)*1000; 

ce1 = reshape(out1.ce, [1 nx*nt]); 
ce2 = reshape(out2.ce, [1 nx*nt]); 
NRMSE.ce = NRMSE_fcn(ce1,ce2)*1000;

cs_bar1 = reshape([out1.cs_bar(1:nn,:); out1.cs_bar(nnp+1:end,:)], [1 nnp*nt]);  
cs_bar2 = reshape([out2.cs_bar(1:nn,:); out2.cs_bar(nnp+1:end,:)], [1 nnp*nt]);  
NRMSE.cs_bar = NRMSE_fcn(cs_bar1,cs_bar2)*1000; 

jn1 = reshape([out1.jn(1:nn,:); out1.jn(nnp+1:end,:)], [1 nnp*nt]);  
jn2 = reshape([out2.jn(1:nn,:); out2.jn(nnp+1:end,:)], [1 nnp*nt]);  
NRMSE.jn = NRMSE_fcn(jn1,jn2)*1000; 

NRMSE.average_internal = mean([NRMSE.phie; NRMSE.ce; NRMSE.cs_bar; NRMSE.jn]); 
end
