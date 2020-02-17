%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example script on how to use the BEST toolbox
%
% This file is a part of the TOOlbox for FAst Battery simulation (TOOFAB)
% Github: https://github.com/Zuan-Khalik/TOOFAB
%
% Author: Zuan Khalik (z.khalik@tue.nl)
%
% TOBASIM is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear all; close all 

%% Prepare parameters
%load current profile
load i_app_WLTP.mat

%Select which parameter set to simulate with. Choice between 'HP' (High
%Power) or 'HE' (High Energy)
parameter_set = 'HE'; 

%Define initial state-of-charge (SOC)
soc_init = 0.95; 
% Define time vector
t = 1:3600; 

if strcmp(parameter_set,'HE')
    Cap = 29.5; 
    p = parameters_LS(); %get parameter struct
else
    Cap = 6; 
    p = parameters_KS(); %get parameter struct
end

% Define input_current in the required format
input_current = [t' 6*(Cap/12.7)*i_app']; 
%% Simulate model
out= DFN(input_current,3600,soc_init,p);
%% Plot results
% Choose at which time step to show internal state variables
eval_time = 1800; 

% Plot the output and some internal state variables
make_plots({out},eval_time,1,1,0)
set(gcf, 'Position',  [20, 20, 800, 950])
set(findall(gcf,'-property','FontSize'),'FontSize',18)
lgd.FontSize = 16; 
