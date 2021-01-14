%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to simulate and plot the considered current profile shown in 
% Fig. 3 in [1] for various selected models as presented in [1] for the 
% high-power (HP) cell [2]. Results shown in Fig. 8 in [1]. 
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
% [2] Smith et al., Control oriented 1d electrochemical model of lithium 
% ion battery, Energy Conversion Management, 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
addpath('Functions')
clear all; close all 
load i_app_validation.mat
i_app = [1.5 0 i_app_validation(1201:end)]; 

grid_high = [10 5 10 20 20]; 
grid_medium = [5 5 5 16 14]; 
grid_low = [3 2 3 9 9];  

soc_init = 1e-6; 
Cap = 7.2; 
% Define time vector
t_interp = [1 1251 1601:4000]; 
i_app_interp = i_app*Cap;

input_current = [t_interp' i_app_interp']; 

N_iter = 11;
%%
p = parameters_KS(grid_high);  
p.current_interp = 'previous';
p.current_extrap = 'nearest'; 
p.set_simp = [1 1 1 1 0 0];
for k = 1:N_iter
    out_CDFN= DFN(input_current,4000,soc_init,p);
    sim_time(k) = out_CDFN.sim_time; 
end
if N_iter>1
    sim_times.CDFN = mean(sim_time(2:end)); 
end

%%
p = parameters_KS(grid_medium); 
p.current_interp = 'previous';
p.current_extrap = 'nearest'; 
p.set_simp = [2 2 2 1 1 0];
for k = 1:N_iter
    out_SDFN_HIFI= DFN(input_current,4000,soc_init,p);
    sim_time(k) = out_SDFN_HIFI.sim_time; 
end
if N_iter > 1
    sim_times.SDFN_HIFI = mean(sim_time(2:end)); 
end

%%
p = parameters_KS(grid_low);  
p.current_interp = 'previous';
p.current_extrap = 'nearest'; 
p.set_simp = [2 2 2 2 1 0];
for k = 1:N_iter
    out_SDFN_LOFI= DFN(input_current,4000,soc_init,p);
    sim_time(k) = out_SDFN_LOFI.sim_time;
end
if N_iter > 1
    sim_times.SDFN_LOFI = mean(sim_time(2:end)); 
end

%%
p = parameters_KS(grid_medium); 
p.current_interp = 'previous';
p.current_extrap = 'nearest'; 
p.set_simp = [2 2 2 1 1 0];
for k = 1:N_iter
    out_SPM_HIFI= SPM(input_current,4000,soc_init,p);
    sim_time(k) = out_SPM_HIFI.solution_time;
end
if N_iter > 1
    sim_times.SPM_HIFI = mean(sim_time(2:end)); 
end

%%
p = parameters_KS(grid_low);  
p.current_interp = 'previous';
p.current_extrap = 'nearest'; 
p.set_simp = [2 2 2 2 1 0];
for k = 1:N_iter
    out_SPM_LOFI= SPM(input_current,4000,soc_init,p);
    sim_time(k) = out_SPM_LOFI.solution_time;
end
if N_iter >1
    sim_times.SPM_LOFI = mean(sim_time(2:end)); 
end
%%
p = parameters_KS([10 10 10 10 10]);  
p.current_interp = 'previous';
p.current_extrap = 'nearest'; 
p.set_simp = [1 1 2 2 0 0];
for k = 1:N_iter
    out_DFN= DFN(input_current,4000,soc_init,p);
    sim_time(k) = out_DFN.sim_time;
end
if N_iter >1
    sim_times.DFN = mean(sim_time(2:end)); 
end
%%
close all
figure(1)
eval_time = 3167; %Show plots at t = 1300 s
make_plots({out_CDFN,out_SDFN_HIFI,out_SDFN_LOFI,out_SPM_HIFI,out_SPM_LOFI},eval_time,1,1,1)
lgd = legend({'CDFN','SDFN-HIFI','SDFN-LOFI','SPM-HIFI','SPM-LOFI'},'Interpreter','latex','FontWeight','bold'); 

set(gcf, 'Position',  [20, 20, 800, 750])
set(findall(gcf,'-property','FontSize'),'FontSize',20)
lgd.FontSize = 16; 
%%
NRMSE.SDFN_HIFI = NRMSE_fcn(out_CDFN.V,out_SDFN_HIFI.V)*1000; 
NRMSE.SDFN_LOFI = NRMSE_fcn(out_CDFN.V,out_SDFN_LOFI.V)*1000; 
NRMSE.SPM_HIFI = NRMSE_fcn(out_CDFN.V,out_SPM_HIFI.V)*1000; 
NRMSE.SPM_LOFI = NRMSE_fcn(out_CDFN.V,out_SPM_LOFI.V)*1000; 
NRMSE.DFN = NRMSE_fcn(out_CDFN.V,out_DFN.V)*1000; 