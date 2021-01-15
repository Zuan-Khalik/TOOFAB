%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to simulate the considered current profile shown in Fig. 3 in [1]
% with the numerical method presented in [2]
%
% Model Simplifications and Their Impact on Computational Complexity for an 
% Electrochemistry-Based Battery Modeling Toolbox
%
% Authors: Z. Khalik, M.C.F. Donkers, H.J. Bergveld
%
% This file is licensed under the BSD 3-Clause License
%
% References
% [1] Khalik et al., Model Simplifications and Their Impact on Computational 
% Complexity for an Electrochemistry-Based Battery Modeling Toolbox, 
% Journal of Power Sources, 2021
% [2] Xia et al, A computationally efficient implementation of a full and
% reduced-order electrochemistry-based model for Li-ion batters, Applied
% Energy, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

addpath('Functions')
clear all; close all 
grid_param = 10*ones(1,5); 
p = parameters_LS(grid_param); 
p.iter_max = 1e4;
load i_app_validation.mat
i_app = [1.5 0 i_app_validation(1201:end)]; 
p.current_interp = 'previous';
p.current_extrap = 'nearest'; 
% grid_param= [10,10,10,10,10];
soc_init = 0; 
Cap = 29.5; 
% Define time vector
t_interp = [1 1251 1601:4000]; 
i_app_interp = i_app*Cap;

input_current = [t_interp' i_app_interp']; 

p.tol = 9e-3;
p.dlnfdx = 0; 
p.set_simp = [0 2 2 2 0 0]; 
p.dt = 1;
%%
tf = 4000;
p.Cap0 = p.Cbat;
p.T_enable = 0;
total_time = tic(); 
out_LX= DFN_LX(input_current,tf,soc_init,p);
out_LX.sim_time_total = toc(total_time); 
sim_time = out_LX.sim_time; 
sim_time_total = out_LX.sim_time_total;

save('Data/results_LX.mat','out_LX')
