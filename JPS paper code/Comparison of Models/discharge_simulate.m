addpath('Models')
clear all; close all 
load i_app_validation.mat
Crate = 0.1; 
i_app = -Crate*ones(1,2); 
% grid_param= [10,10,10,10,10];
grid_param = 10*ones(1,5); 
soc_init = 1; 
Cap = 29.5; 
% Define time vector
t_interp = [1 1e6]; 
i_app_interp = i_app*Cap;

input_current = [t_interp' i_app_interp']; 

N_iter = 1; 
%%
p = parameters_LS([10 8 13 3 3]);  
p.set_simp = [2 2 2 2 0 0];
p.dt = p.dt/Crate; 
out_CDFN= DFN(input_current,1e6,soc_init,p);

%%
p = parameters_LS([10 8 13 3 3]); 
p.set_simp = [2 2 2 1 0 0];
for k = 1:N_iter
    out_SDFN_HIFI= DFN(input_current,4000,soc_init,p);
    sim_time(k) = out_SDFN_HIFI.solution_time; 
end
if N_iter > 1
    simtime_SDFN_HIFI = mean(sim_time(2:end)); 
end
%%
p = parameters_LS([6 2 8 3 3]);  
p.set_simp = [2 2 2 2 0 0];
for k = 1:N_iter
    out_SDFN_LOFI= DFN(input_current,4000,soc_init,p);
    sim_time(k) = out_SDFN_LOFI.solution_time;
end
if N_iter > 1
    simtime_SDFN_LOFI = mean(sim_time(2:end)); 
end
%%
p = parameters_LS([10 8 13 3 3]); 
p.set_simp = [2 2 2 1 0 0];
for k = 1:N_iter
    out_SPM_HIFI= SPM(input_current,4000,soc_init,p);
    sim_time(k) = out_SPM_HIFI.solution_time;
end
if N_iter > 1
    simtime_SPM_HIFI = mean(sim_time(2:end)); 
end
%%
p = parameters_LS([6 2 8 3 3]);  
p.set_simp = [2 2 2 2 0 0];
for k = 1:N_iter
    out_SPM_LOFI= SPM(input_current,4000,soc_init,p);
    sim_time(k) = out_SPM_LOFI.solution_time;
end
if N_iter >1
    simtime_SPM_LOFI = mean(sim_time(2:end)); 
end
%%
close all
figure(1)
eval_time = 3167; %Show plots at t = 1300 s
make_plots({out_CDFN,out_SDFN_HIFI,out_SDFN_LOFI,out_SPM_HIFI,out_SPM_LOFI},eval_time,1,1,0)

set(gcf, 'Position',  [20, 20, 1000, 750])
set(findall(gcf,'-property','FontSize'),'FontSize',16)

%%
% NRMSE_DFN = NRMSE_fcn(out_CDFN.V,out_DFN.V); 
NRMSE_SDFN_HIFI = NRMSE_fcn(out_CDFN.V,out_SDFN_HIFI.V)*1000; 
NRMSE_SDFN_LOFI = NRMSE_fcn(out_CDFN.V,out_SDFN_LOFI.V)*1000; 
NRMSE_SPM_HIFI = NRMSE_fcn(out_CDFN.V,out_SPM_HIFI.V)*1000; 
NRMSE_SPM_LOFI = NRMSE_fcn(out_CDFN.V,out_SPM_LOFI.V)*1000; 
