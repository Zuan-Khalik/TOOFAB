%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to simulate Fig. 5 in [1] for the high-power (HP) cell [2]. 
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
% [2] Smith et al., Control oriented 1d electrochemical model of lithium 
% ion battery, Energy Conversion Management, 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
addpath('Functions')
clear all; close all;

Crate = 2;
soc_init = 0.5; 
dt = 1;
time = 1200;
t = dt:dt:time;
load('i_app.mat')
i_app = 1.5*(1/6)*7.2*i_app(1:1200); 
input_current = [t' i_app']; 
p.set_simp = [2 2 2 2 1 0]; 
p.options.verbose = 1; 

grid_base = 40*ones(1,5);
p = parameters_KS(grid_base); 
n_iter = 11; 
for i = 1:n_iter
    out = DFN(input_current,time,soc_init,p);
    V_base = out.V;
    ce_base = out.ce;
    cs_base = out.stoich;
    jn_base = out.jn; 
    x_base = out.x;
    sim_time_base_i(i) = out.sim_time; 
end
sim_time_base = mean(sim_time_base_i(2:end)); 
% % 

%%
    grid_select{1} = [3 2 3 9 9]; 
    grid_select{2} = [5 5 5 16 14]; 
    grid_select{3} = [10 5 10 20 20];
    
for j = 1:length(grid_select)
    model_order(j) = grid_select{j}(1)+grid_select{j}(2)+grid_select{j}(3)+grid_select{j}(1)*grid_select{j}(4)+grid_select{j}(3)*grid_select{j}(5); 
    for k = 1:n_iter
        p = parameters_KS(grid_select{j}); 
        out = DFN(input_current,time,soc_init,p);
        V{j} = out.V;
        ce{j} = out.ce; 
        cs{j} = out.stoich;
        jn{j} = out.jn; 
        x{j} = out.x; 
        sim_time(j,k) = out.sim_time;
        rmserror(j) = NRMSE_fcn(V_base,V{j}); 
    end
    sim_time_mean(j) = mean(sim_time(j,2:end),2)'; 
end

save('Data/compare_model_orders_HP.mat','x','x_base','cs_base','cs','ce_base','ce','jn_base','jn','V_base','V','sim_time_mean','rmserror','t','grid_select')