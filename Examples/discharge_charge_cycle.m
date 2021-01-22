%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example showing how the toolbox can be used to control the battery in
% closed-loop, where the battery is simulated with a discharge-charge cycle
%
% This file is a part of the TOOlbox for FAst Battery simulation (TOOFAB)
% Github: https://github.com/Zuan-Khalik/TOOFAB
%
% Author: Zuan Khalik (z.khalik@tue.nl)
%
% TOOFAB is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
addpath('Functions')

clear all; close all; 

%Apply Simplification [S1] and [S2-II]-all
p.set_simp = [2 2 2 2 1 0];

%Choose C-rate for charging and discharging
p.Crate_charge = 0.5; 
p.Crate_discharge = 0.5;

%Choose step sizes for the various phases in the cycle. Note that the
%parameter values have no significant meaning, and is only meant to show
%how the additional dt output of the input_current function can be used
p.dt_charge = 4;
p.dt_discharge = 4;
p.dt_CV = 3;

%Run simulation 
out = DFN(@cycle,2e5,1,p); 

%% Plot results
% Choose at which time step to show internal state variables
eval_time = 1800; 

% Plot the output and some internal state variables
make_plots({out},eval_time,1,1,0)
set(gcf, 'Position',  [20, 20, 800, 950])
set(findall(gcf,'-property','FontSize'),'FontSize',18)
lgd.FontSize = 16; 


function [i_app_next,mem,end_simulation,dt] = cycle(k,t,i_app,V,soc,cs,ce,phis,phie,mem,p)
Cbat = 1*p.Cbat/3600; %Capacity of the battery
V_max =p.OCV_max;   %Voltage at which the battery switches to CV mode
V_min =p.OCV_min;   %Voltage after which the battery stops discharging
i_min = 0.05*Cbat;  %Current level to terminate the CV mode with
Kp1 = 50;   %Proportional gain of the PI controller for the CV stage
Ki1 = 10;   %Integral gain for the PI controller

end_simulation=0; %flag to indicate whether the simulation should end after the iteration. Has to be set at every iteration

% Initialize the various parameters. Note how the mem struct is used to
% store the states of the PI controller 
if k==1
    mem.status = 0;
    mem.e = 0;
    mem.e2 = 0;
end

% Condition to terminate the initial CC discharge stage and start the CC
% charge stage
if V(k)<=V_min && mem.status==0
    mem.status = 1;  
end

% Condition to terminate the CC charge stage and start the CV stage
if V(k)>=V_max &&mem.status==1
    mem.status = 2; 
end

%Condition to terminate the simulation
if i_app(k)<i_min && mem.status==2
    end_simulation=1; 
end

%Switch between the different stages
switch mem.status
    case 0 %CC discharge
        i_app_next = -p.Crate_discharge*Cbat; %value of i_app for the next iteration
        dt = p.dt_discharge; %assign the chosen step size
    case 1 %CC charge
        i_app_next = p.Crate_charge*Cbat; %value of i_app for the next iteration
        dt = p.dt_charge; %assign the chosen step size
    case 2 %CV
        %Solve the PI controller system (in this case with a forward Euler
        %discretization
        mem.e(k+1) = V_max-V(k); 
        mem.e2(k+1) = (mem.e(k+1)-mem.e(k))/p.dt; 
        i_app_next = i_app(k)+p.dt*(Kp1*mem.e2(k+1)+Ki1*mem.e(k+1));
        dt = p.dt_CV; %assign the chosen step size
end
end