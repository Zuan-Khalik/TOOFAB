%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Model Simplifications and Their Impact on Computational Complexity for an 
% Electrochemistry-Based Battery Modeling Toolbox
%
% Authors: Z. Khalik, M.C.F. Donkers, H.J. Bergveld
%
% This file is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
addpath('Functions')
clear all; close all 
grid_param= 10*ones(1,5); 
p = parameters_KS(grid_param); 
load i_app_validation.mat
i_app = [1.5 0 i_app_validation(1201:end)]; 
soc_init = 1e-6;

Cap = 7.2; 
t_interp = [1 1251 1601:4000]; 
i_app_interp = i_app*Cap;
p.current_interp = 'previous';
p.current_extrap = 'nearest'; 

input_current = [t_interp' i_app_interp']; 

n_iter = 11;

set_simp.full = [0 0 0 0 0 0]; 
set_simp.s1 =[0 0 0 0 1 0]; 
set_simp.s22_kappa = [2 0 0 0 0 0]; 
set_simp.s22_De = [0 2 0 0 0 0]; 
set_simp.s22_nu = [0 0 2 0 0 0]; 
set_simp.s22_Ds = [0 0 0 2 0 0];
set_simp.s22 = [2 2 2 2 0 0]; 
set_simp.s21 = [1 1 1 1 0 0]; 

names = {'full', 's1','s22_kappa', 's22_De', 's22_nu', 's22_Ds', 's22', 's21'}; 

for i = 1:length(names)
    p.set_simp = set_simp.(names{i}); 
    for k = 1:n_iter
        out.(names{i}) = DFN(input_current,4000,soc_init,p);
        simtime.(names{i})(k) = out.(names{i}).solution_time; 
    end
end

%%
for i = 1:length(names)
    RMSE.(names{i}) =  compute_nrmse(out.(names{i}),out.full); 
    simtime_mean.(names{i}) = mean(simtime.(names{i})(2:end)); 
end

save('Data/simplifications_data_HP.mat','RMSE','simtime_mean','names')