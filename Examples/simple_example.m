%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example script on how to use the BEST toolbox
%
% This file is a part of the TOOlbox for BAttery SIMulation (TOBASIM)
% Github: https://github.com/Zuan-Khalik/Battery-Simulation-Toolbox
%
% Author: Zuan Khalik (z.khalik@tue.nl)
%
% TOBASIM is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear all; close all 
addpath('Functions')
%% Prepare parameters
%Select which parameter set to simulate with. Choice between 'HP' (High
%Power) or 'HE' (High Energy)
parameter_set = 'HE'; 

%Define initial state-of-charge (SOC)
soc_init = 0; 
% Define time vector
t = 1:3600; 

if strcmp(parameter_set,'HE')
    Cap = 29.5; 
    p = parameters_LS(); %get parameter struct
else
    Cap = 6; 
    p = parameters_KS(); %get parameter struct
end
p.verbose = 2;
% Define input_current in the required format
load i_app_validation.mat
i_app = [1.5 1.5 0 0 i_app_validation(1201:end)]; 
t_interp = [1 1250 1251 1600 1601:4000]; 
i_app_interp = i_app*Cap;

input_current = [t_interp' i_app_interp']; 

%% Simulate model
out= DFN(input_current,4000,soc_init,p);

%% Plot results
% Choose at which time step to show internal state variables
eval_time = 1800; 

% Plot the output and some internal state variables
make_plots({out},eval_time,1,1,0)
set(gcf, 'Position',  [20, 20, 800, 950])
set(findall(gcf,'-property','FontSize'),'FontSize',18)
lgd.FontSize = 16; 
