%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to convert the LIONSIMBA [2] output file to the form used in the
% TOOFAB toolbox [1] for plotting
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
clear all; close all 
load i_app_validation.mat
i_app = [1.5 0 i_app_validation(1201:end)]; 
grid_param = 10*ones(1,5); 
soc_init = 0; 
Cap = 29.5; 
% Define time vector
t_interp = [1 1251 1601:4000]; 
i_app_interp = i_app*Cap;

input_current = [t_interp' i_app_interp']; 

p = parameters_LS(grid_param); 
p.set_simp = [1 2 2 2 0 0]; 
p.dt = 1;
%%
total_time = tic(); 
out_ZK= DFN(input_current,4000,soc_init,p);
out_ZK.total_time = toc(total_time); 
%%
load('Data/results_LS_orig.mat')
%%
out_LS.V = out.Voltage{1}; 
out_LS.t = out.time{1};
out_LS.x = out_ZK.x; 
out_LS.p = out_ZK.p;
out_LS.sim_time = out.simulation_time; 
out_LS.sim_time_total = out.sim_time_total; 
n_t = length(out_LS.V); 
out_LS.phie = flip(out.Phie{1}-out.Phie{1}(:,1),2); 
out_LS.ce = flip(out.ce{1},2); 
out_LS.U = flip([out.Up{1} NaN(n_t,p.ns) out.Un{1}],2); 
out_LS.phis = flip([out.Phis{1}(:,1:p.nn) NaN(n_t,p.ns) out.Phis{1}(:,p.nn+1:end)]-out.Phie{1}(:,1),2);
out_LS.cs_bar = flip([out.cs_surface{1}(:,1:p.nn) NaN(n_t,p.ns) out.cs_surface{1}(:,out_ZK.p.nn+1:end)],2);
out_LS.stoich = flip([out.cs_surface{1}(:,1:p.nn)/out_ZK.p.cs_max_pos NaN(n_t,p.ns) out.cs_surface{1}(:,out_ZK.p.nn+1:end)/out_ZK.p.cs_max_neg],2);
out_LS.eta = flip([out.etap{1} NaN(n_t,p.ns) out.etan{1}],2);
p.k0 = out_ZK.p.k0; 
p.cs_bar_max = out_ZK.p.cs_bar_max; 
elec_range = out_ZK.p.elec_range; 
cs = out_LS.cs_bar'; 
ce = out_LS.ce'; 
eta = out_LS.eta';
T = out_ZK.T'; 
for k = 1:4000
    deltap = (0.5*p.F)./(p.R*T(k)).*eta(elec_range,k);
    ip = 2*p.k0.*sqrt(ce(elec_range,k)).*sqrt(cs(elec_range,k)).*sqrt(p.cs_bar_max-cs(elec_range,k));
    jn(:,k) = (1/p.F)*ip.* sinh(deltap);
end
out_LS.jn = [jn(1:p.nn,:); NaN(p.ns,n_t); jn(p.nn+1:end,:)]';
save('results_LS.mat')
