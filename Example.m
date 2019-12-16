%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example script on how to use the BEST toolbox
%
% This file is a part of the BattEry Simulation Toolbox (BEST)
% Github: https://github.com/Zuan-Khalik/Battery-Simulation-Toolbox
%
% Author: Zuan Khalik (z.khalik@tue.nl)
%
% BEST is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear all; close all 

%% Prepare parameters
%load current profile
load i_app_WLTP.mat

%Select which parameter set to simulate with. Choice between 'HP' (High
%Power) or 'HE' (High Energy)
parameter_set = 'HP'; 

%Define initial state-of-charge (SOC)
soc_init = 0.8; 
% Define time vector
t = 1:3600; 

if strcmp(parameter_set,'HE')
    Cap = 12.7; 
    p = parameters_JN(); %get parameter struct
else
    Cap = 6; 
    p = parameters_KS(); %get parameter struct
end

% Define input_current in the required format
input_current = [t' (12.7/Cap)*i_app']; 

%% Simulate model
out= DFN(input_current,3600,soc_init,p);

%% Plot results
% Choose at which time step to show internal state variables
eval_time = 1800; 

% Plot the output and some internal state variables
make_plots(out,eval_time,'k','k',1,0)
