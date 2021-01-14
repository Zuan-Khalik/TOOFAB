%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to simulate the considered current profile shown in Fig. 3 in [1]
% with the numerical method presented in [1]
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
addpath('Functions')
clear all; close all 
grid_param = 10*ones(1,5); 
p = parameters_LS(grid_param); 
load i_app_validation.mat
i_app = [1.5 1.5 0 0 i_app_validation(1201:end)]; 
p.current_interp = 'previous';
p.current_extrap = 'nearest'; 

soc_init = 0; 
Cap = 29.5; 
% Define time vector
t_interp = [1 1250 1251 1600 1601:4000]; 
i_app_interp = i_app*Cap;

input_current = [t_interp' i_app_interp']; 

p.tol = 3e-2;
% p.tol = 1e-6;
p.dlnfdce = 0; 
p.dlnfdx = 0;
p.set_simp = [0 2 2 2 0 0]; 
p.dt = 1;
%%

tf = 4000;

total_time = tic(); 
out_ZK= DFN(input_current,tf,soc_init,p);
out_ZK.sim_time_total = toc(total_time); 
sim_time = out_ZK.sim_time; 
sim_time_total =out_ZK.sim_time_total;

save('Data/results_ZK.mat','out_ZK')