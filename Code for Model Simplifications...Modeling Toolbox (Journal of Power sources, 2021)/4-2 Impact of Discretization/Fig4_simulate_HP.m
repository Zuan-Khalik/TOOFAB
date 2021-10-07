%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to simulate Fig. 4 in [1] for the high-power (HP) cell [2]. 
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
clear all; close all;
addpath('Functions')
Crate = 2;
soc_init = 0.5; 
dt = 1;
time = 1200;
t = dt:dt:time;
load('i_app.mat')
i_app = 1.5*(1/6)*7.2*i_app(1:1200); 
input_current = [t' i_app']; 
p = parameters_KS();
p.set_simp = [2 2 2 2 1 0]; 
p.verbose = 0; 
n_rep = 11; % number of repetitions
%------- Uncomment below to generate data_grid_test dataset--------------
% load data_grid_test.mat

grid_select = [38,34,28,24,20,18,16, 14,12,10,9,8,7,6, 5,4,3,2]; 
grid_base = {40, 40 ,40 , 40, 40};

for i = 1:n_rep
    [p.grid.nn,p.grid.ns,p.grid.np,p.grid.nrp,p.grid.nrn] = grid_base{:}; 
    out = DFN(input_current,time,soc_init,p);
    V_base_temp{i} = out.V;
    sim_time_base_i(i) = out.sim_time; 
end
V_base = V_base_temp{1}; 
sim_time_base = mean(sim_time_base_i); 
% % 
for i = 1:5
    for j = 1:length(grid_select)
        for k = 1:n_rep
            grid_loop_temp = grid_base;  
            grid_loop_temp{i} = grid_select(j); 
            grid_loop{i,j,k} = grid_loop_temp; 
        end
    end
end

f = waitbar(0);
for i = 1:5
    for j = 1:length(grid_select)
        waitbar(((i-1)*length(grid_select)+j)/(5*length(grid_select)),f,['i = ', num2str(i),'/5, ','j = ',num2str(j),'/',num2str(length(grid_select))])
        for k = 1:n_rep
            [p.grid.nn,p.grid.ns,p.grid.np,p.grid.nrp,p.grid.nrn] = grid_loop{i,j,k}{:}; 
            out = DFN(input_current,time,soc_init,p);
            V_temp2{k} = out.V; 
            sim_time_temp2(k) = out.sim_time;
        end
        V_temp{j} = V_temp2{1}; 
        sim_time_temp(j,:) = sim_time_temp2; 
        rmserror_temp(j) = NRMSE_fcn(V_base,V_temp{j});
    end
    sim_time{i} = sim_time_temp; 
    V{i,:} = V_temp;
    rmserror{i} = rmserror_temp; 
end

for i = 1:5
    rmserror_mean{i} = mean(rmserror{i},2)'; 
    sim_time_mean{i} = mean(sim_time{i},2)'; 
end

%%
figure(1)
subplot(2,1,1)
semilogy([grid_base{1}, grid_select],[0,rmserror{1}]*1000,'k','LineWidth',2)
hold on
semilogy([grid_base{1}, grid_select],[0,rmserror{2}]*1000,'k:','LineWidth',2)
hold on
semilogy([grid_base{1}, grid_select],[0,rmserror{3}]*1000,'k--','LineWidth',2)
hold on
semilogy([grid_base{1}, grid_select],[0,rmserror{4}]*1000,'r','LineWidth',2)
hold on
semilogy([grid_base{1}, grid_select],[0,rmserror{5}]*1000,'r--','LineWidth',2)
xlabel('Value of varying grid parameter')
ylabel('RMS error [mV]')
grid on

subplot(2,1,2)
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{1}],'k','LineWidth',2)
hold on
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{2}],'k:','LineWidth',2)
hold on
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{3}],'k--','LineWidth',2)
hold on
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{4}],'r','LineWidth',2)
hold on
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{5}],'r--','LineWidth',2)
xlabel('Value of varying grid parameter')
ylabel('Simulation time [s]')
grid on

legend('Vary n_n', 'Vary n_s', 'Vary n_p', 'Vary n_{r,n}', 'Vary n_{r,p}')

%%
save('Data/impact_discretization_HP2','t','i_app','grid_base','grid_select','rmserror','sim_time_mean','sim_time_base','soc_init','time','V','V_base')