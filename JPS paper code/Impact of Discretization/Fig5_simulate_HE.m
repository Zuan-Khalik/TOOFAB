%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to simulate Fig. 5 in [1] for the high-energy (HE) cell [2]. 
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
% [2] Torchio et al., A matlab framework based on a finite novolume model
% suitable for Li-ion battery design, simulation, and control, Journal of
% the Electrochemical Society, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
addpath('Functions')
clear all; close all;

Crate = 2;
soc_init = 0.5; 
dt = 1;
time = 1200;
t = dt:dt:time;
load('i_app.mat')
i_app = 1.5*(1/6)*29.5*i_app(1:1200); 
input_current = [t' i_app']; 
p.set_simp = [2 2 2 2 1 0]; 
p.options.verbose = 1; 
%------- Uncomment below to generate data_grid_test dataset--------------
% load data_grid_test.mat

grid_base = 40*ones(1,5);
p = parameters_LS(grid_base); 
n_iter = 4; 
for i = 1:n_iter
    out = DFN(input_current,time,soc_init,p);
    V_base = out.V;
    ce_base = out.ce;
    cs_base = out.stoich;
    jn_base = out.jn; 
    x_base = out.x;
    sim_time_base_i(i) = out.solution_time; 
end
sim_time_base = mean(sim_time_base_i(2:end)); 
% % 

%%
    grid_select{1} = [4,2,6,3,3]; 
    grid_select{2} = [12,8,15,3,3]; 
    grid_select{3} = [22 10 23 3 7];
    
for j = 1:length(grid_select)
    model_order(j) = grid_select{j}(1)+grid_select{j}(2)+grid_select{j}(3)+grid_select{j}(1)*grid_select{j}(4)+grid_select{j}(3)*grid_select{j}(5); 
    for k = 1:n_iter
        p = parameters_LS(grid_select{j}); 
        out = DFN(input_current,time,soc_init,p);
        V{j} = out.V;
        ce{j} = out.ce; 
        cs{j} = out.stoich;
        jn{j} = out.jn; 
        x{j} = out.x; 
        sim_time(j,k) = out.solution_time;
        rmserror(j) = NRMSE_fcn(V_base,V{j}); 
    end
    sim_time_mean(j) = mean(sim_time(j,2:end),2)'; 
end

save('Data/compare_model_orders_HE.mat','x','x_base','cs_base','cs','ce_base','ce','jn_base','jn','V_base','V','sim_time_mean','rmserror','t','grid_select')