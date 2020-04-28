clear all; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to obtain the data for the  comparison of the numerical methods 
% of [1],[2],[3] Fig. 6 in the associated paper [1]) 
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
% [2] Xia et al, A computationally efficient implementation of a full and
% reduced-order electrochemistry-based model for Li-ion batters, Applied
% Energy, 2017
% [3] Torchio et al., A matlab framework based on a finite novolume model
% suitable for Li-ion battery design, simulation, and control, Journal of
% the Electrochemical Society, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Obtain LIONSIMBA data 
% Please make sure that LIONSIMBA is correctly installed
LIONSIMBA_simulation
restructure_LIONSIMBA

%% Obtain LX data
LU_simulation

%% Obtain ZK data
ZK_simulation