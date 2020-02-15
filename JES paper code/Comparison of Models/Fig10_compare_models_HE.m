%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to simulate and plot the considered current profile shown in 
% Fig. 4 in [1] for various selected models as presented in [1] for the 
% high-energy (HE) cell [2]. Results shown in Fig. 10 in [1]. 
%
% On Trade-offs between Computational Complexity and Accuracy of
% Electrochemistry-based Battery models
%
% Authors: Z. Khalik, H.J. Bergveld, M.C.F. Donkers
%
% This file is licensed under the BSD 3-Clause License
%
% References
% [1] Khalik et al., On trade-offs between Computational Complexity and 
% Accuracy of Electrochemistry-based Battery Models, Journal of the 
% Electrochemical Society, 2020, submitted
% [2] Torchio et al., A matlab framework based on a finite novolume model
% suitable for Li-ion battery design, simulation, and control, Journal of
% the Electrochemical Society, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
addpath('Functions')
clear all; close all 
load i_app_validation.mat
i_app = [1.5 0 i_app_validation(1201:end)]; 

grid_high = [22 10 23 3 7]; 
grid_medium = [12 8 15 3 3]; 
grid_low = [9 2 12 3 3];  

soc_init = 0; 
Cap = 29.5; 
% Define time vector
t_interp = [1 1251 1601:4000]; 
i_app_interp = i_app*Cap;

input_current = [t_interp' i_app_interp']; 

N_iter = 4; 
%%
p = parameters_LS(grid_high);  
p.set_simp = [1 1 1 1 0 0];
for k = 1:N_iter
    out_CDFN= DFN(input_current,4000,soc_init,p);
    sim_time(k) = out_CDFN.solution_time; 
end
if N_iter>1
    simtime_CDFN = mean(sim_time(2:end)); 
end
%%
p = parameters_LS(grid_medium); 
p.set_simp = [2 2 2 1 0 1];
for k = 1:N_iter
    out_SDFN_HIFI= DFN(input_current,4000,soc_init,p);
    sim_time(k) = out_SDFN_HIFI.solution_time; 
end
if N_iter > 1
    simtime_SDFN_HIFI = mean(sim_time(2:end)); 
end
%%
p = parameters_LS(grid_low);  
p.set_simp = [2 2 2 2 1 0];
for k = 1:N_iter
    out_SDFN_LOFI= DFN(input_current,4000,soc_init,p);
    sim_time(k) = out_SDFN_LOFI.solution_time;
end
if N_iter > 1
    simtime_SDFN_LOFI = mean(sim_time(2:end)); 
end
%%
p = parameters_LS(grid_medium); 
p.set_simp = [2 2 2 1 0 0];
for k = 1:N_iter
    out_SPM_HIFI= SPM(input_current,4000,soc_init,p);
    sim_time(k) = out_SPM_HIFI.solution_time;
end
if N_iter > 1
    simtime_SPM_HIFI = mean(sim_time(2:end)); 
end
%%
p = parameters_LS(grid_low);  
p.set_simp = [2 2 2 2 1 0];
for k = 1:N_iter
    out_SPM_LOFI= SPM(input_current,4000,soc_init,p);
    sim_time(k) = out_SPM_LOFI.solution_time;
end
if N_iter >1
    simtime_SPM_LOFI = mean(sim_time(2:end)); 
end
%%
p = parameters_LS([10 10 10 10 10]);  
p.set_simp = [1 1 2 2 0 0];
for k = 1:N_iter
    out_DFN= DFN(input_current,4000,soc_init,p);
    sim_time(k) = out_DFN.solution_time;
end
if N_iter >1
    simtime_DFN = mean(sim_time(2:end)); 
end
%%
close all
figure(1)
eval_time = 3167; %Show plots at t = 1300 s
make_plots({out_CDFN,out_SDFN_HIFI,out_SDFN_LOFI,out_SPM_HIFI,out_SPM_LOFI},eval_time,1,1,0)
lgd = legend({'CDFN','SDFN-HIFI','SDFN-LOFI','SPM-HIFI','SPM-LOFI'},'Interpreter','latex','FontWeight','bold'); 

set(gcf, 'Position',  [20, 20, 800, 750])
set(findall(gcf,'-property','FontSize'),'FontSize',20)
lgd.FontSize = 16; 
%%
NRMSE_SDFN_HIFI = NRMSE_fcn(out_CDFN.V,out_SDFN_HIFI.V)*1000; 
NRMSE_SDFN_LOFI = NRMSE_fcn(out_CDFN.V,out_SDFN_LOFI.V)*1000; 
NRMSE_SPM_HIFI = NRMSE_fcn(out_CDFN.V,out_SPM_HIFI.V)*1000; 
NRMSE_SPM_LOFI = NRMSE_fcn(out_CDFN.V,out_SPM_LOFI.V)*1000; 
NRMSE_DFN = NRMSE_fcn(out_CDFN.V,out_DFN.V)*1000;
 