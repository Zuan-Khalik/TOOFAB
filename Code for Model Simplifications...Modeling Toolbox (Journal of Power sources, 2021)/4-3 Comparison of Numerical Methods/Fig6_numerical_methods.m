%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to plot the comparison of the numerical methods of [1],[2],[3]
% (Fig. 6 in the associated paper [1]) 
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
% [3] Torchio et al., A matlab framework based on a finite novolume model
% suitable for Li-ion battery design, simulation, and control, Journal of
% the Electrochemical Society, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

addpath('Functions')
addpath('Data')
clear all; close all

% To obtain the data for plotting this figure, please run
% "Fig8_obtain_data.m" 
load Data/results_LS.mat % Data from simulation using method of [3]
load Data/results_LX.mat % Data from simulation using method of [2]
load Data/results_ZK.mat % Data from simulation using method of [1]
figure(1)
eval_time = 2967; %Show plots at t = eval_time

out_LS.param = out_ZK.param;
make_plots({out_LS,out_LX,out_ZK},eval_time,1,1,0)
lgd = legend({'LIONSIMBA [29]','Method of Xia et al. [13]','Proposed method'},'Interpreter','latex','FontWeight','bold');

set(gcf, 'Position',  [20, 20, 800, 750])
set(findall(gcf,'-property','FontSize'),'FontSize',20)
lgd.FontSize = 16; 

%% NRSME and simulation times
NRMSE.LX = compute_nrmse(out_LS,out_LX);
NRMSE.ZK = compute_nrmse(out_LS,out_ZK);
sim_times.ZK.sim_time_solution = out_ZK.sim_time; 
sim_times.ZK.sim_time_total = out_ZK.sim_time_total; 
sim_times.LS.sim_time_solution = out_LS.sim_time; 
sim_times.LS.sim_time_total = out_LS.sim_time_total; 
sim_times.LX.sim_time_solution = out_LX.sim_time; 
sim_times.LX.sim_time_total = out_LX.sim_time_total; 