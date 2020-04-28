%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to simulate the considered current profile shown in Fig. 3 in [1]
% with the numerical method presented in [2]
%
% Model Simplifications and Its Impact on Computational Complexity for an 
% Electrochemistry-Based Battery Modeling Toolbox
%
% Authors: Z. Khalik, M.C.F. Donkers, H.J. Bergveld
%
% This file is licensed under the BSD 3-Clause License
%
% References
% [1] Khalik et al., Model Simplifications and Its Impact on Computational 
% Complexity for an Electrochemistry-Based Battery Modeling Toolbox, 
% Journal of Power Sources, 2020, submitted
% [2] Xia et al, A computationally efficient implementation of a full and
% reduced-order electrochemistry-based model for Li-ion batters, Applied
% Energy, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

addpath('Functions')
clear all; close all 
load i_app_validation.mat
i_app = [1.5 0 i_app_validation(1201:end)]; 
% grid_param= [10,10,10,10,10];
grid_param = 10*ones(1,5); 
soc_init = 0; 
Cap = 29.5; 
% Define time vector
t_interp = [1 1251 1601:4000]; 
i_app_interp = i_app*Cap;

input_current = [t_interp' i_app_interp']; 

p = parameters_LS(grid_param); 
p.dlnfdx = 0; 
p.set_simp = [0 2 2 2 0 0]; 
p.dt = 1;
%%
total_time = tic(); 
out_LX= DFN_LX(input_current,4000,soc_init,p);
out_LX.sim_time_total = toc(total_time); 
save('results_LX.mat','out_LX')
