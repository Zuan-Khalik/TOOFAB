%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to obtain the data to plot the discharge curves as shown in Fig. 9
% in [1] for the high-power (HP) [3] and high-energy (HE) [2] parameters
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
% [3] Smith et al., Control oriented 1d electrochemical model of lithium 
% ion battery, Energy Conversion Management, 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear all; close all

%% Simulate discharge curves of HP cell 
discharge_simulate_HP

%% Simulate discharge curves of HE cell
discharge_simulate_HE