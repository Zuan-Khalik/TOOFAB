%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to simulate the considered current profile shown in Fig. 3 in [1]
% with the LIONSIMBA toolbox [2]. Please note that some of the code has
% been copied from the LIONSIMBA toolbox. 
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

% LIONSIMBA example script
% Custom current profile: this script shows how to run a custom current
% profile charge/discharge

%   This file is part of the LIONSIMBA Toolbox
%
%	Official web-site: 	http://sisdin.unipv.it/labsisdin/lionsimba.php
% 	Official GitHUB: 	https://github.com/lionsimbatoolbox/LIONSIMBA
%
%   LIONSIMBA: A Matlab framework based on a finite volume model suitable for Li-ion battery design, simulation, and control
%   Copyright (C) 2016-2018 :Marcello Torchio, Lalo Magni, Davide Raimondo,
%                            University of Pavia, 27100, Pavia, Italy
%                            Bhushan Gopaluni, Univ. of British Columbia, 
%                            Vancouver, BC V6T 1Z3, Canada
%                            Richard D. Braatz, 
%                            Massachusetts Institute of Technology, 
%                            Cambridge, Massachusetts 02142, USA
%   
%   Main code contributors to LIONSIMBA 2.0:
%                           Ian Campbell, Krishnakumar Gopalakrishnan,
%                           Imperial college London, London, UK
%
%   LIONSIMBA is a free Matlab-based software distributed with an MIT
%   license.

% Clear the workspace
clear all; close all

%% Create interpolant for current
load i_app_validation.mat
i_app = [1.5 1.5 0 0 i_app_validation(1201:end)]; 
t_interp = [1 1250 1251 1600 1601:4000]; 
i_app_interp = i_app*29.5;
F_current = griddedInterpolant(t_interp,i_app_interp,'linear'); 
save('Data/F_i_app','F_current')

fprintf('Please copy getInputCurrentDensity.m and F_i_app.mat from /Data to /LIONSIMBA-master/battery_model_files/external_functions \n')
fprintf('To install LIONSIMBA, please follow the instructions on the Github page: https://github.com/lionsimbatoolbox/LIONSIMBA \n')
check = input('Press ENTER if this has been done'); 

%% Parameters
% Define the integration times.
t0 = 0;
tf = 10^4;

% Define the initial state structure
initialState.Y  = [];
initialState.YP = [];

% Define the parameters structure.
param{1}               = Parameters_init(0);

param{1}.Np            = 10;
param{1}.Ns            = 10;
param{1}.Nn            = 10;

param{1}.hcell         = 1;
param{1}.Tref          = 298.15;

param{1}.AbsTol        = 1e-6;
param{1}.RelTol        = 1e-6;

param{1}.CutoffSOC     = 0;

param{1}.SolidPhaseDiffusion = 3;
param{1}.TemperatureEnabled = 0;

if(param{1}.SolidPhaseDiffusion == 3)
    multp = param{1}.Nr_p;
    multn = param{1}.Nr_n;
else
    multp = 1;
    multn = 1;
end
I1C = 29.5;

C_rate = 1.5;
%% Discharge section

% Apply a custom current profile
param{1}.OperatingMode = 4;
param{1}.edge_values = 2;

% Run the simulation
sim_time_total_tic = tic;   
out = startSimulation(1,4000,initialState,0,param);
out.sim_time_total = toc(sim_time_total_tic);
out_LS = out; 

save('Data/results_LS','out_LS')