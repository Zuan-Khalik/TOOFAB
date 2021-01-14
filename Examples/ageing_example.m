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
% TOBASIM is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
addpath('Functions')
clear all; close all; 

%Apply Simplification [S1] and [S2-II]-all
p.set_simp = [2 2 2 2 1 0];

%Enable ageing dynamics
p.ageing = 1;

%Disable verbosity 
p.verbose = 0;

%Choose C-rate for charging and discharging
p.charge_current = 30; 
p.discharge_current = 6;

%Choose step sizes for the various phases in the cycle. These step sizes
%have been selected such that a small discretization error is incurred 
p.dt_discharge = 10;
p.dt_charge = 5;
p.dt_CV = 2;

%Initialize loop variables
i = 1;
SoH = 1; %state of health
while SoH(i)>0.8 %terminate when state of health is at 80% 
        %At the first cycle, the initial conditions are static 
        %corresponding to SoC = 1
        if i==1 
            out = DFN(@CC_CV,2e5,1,p); 
            Cbat_fresh = out.param.Cbat;
        %At the remaining cycles we would like to use the states of at the 
        %end of the previous simulation as the initial condition for this simulation    
        else 
            %%
            out = DFN(@CC_CV,2e5,states_init,p); 
        end
        % Save variables
        SoH(i+1) = out.ageing.Cbat_aged/out.param.Cbat; %compute state-of-health
        ageing.Rf{i} = out.ageing.Rf; %SEI film resistance
        ageing.ct(i) = (out.t(out.mem.charge_end)-out.t(out.mem.discharge_end))/60; %charging time 
        ageing.Cbat(i) = out.ageing.Cbat_aged; %Capacity of the aged battery
        ageing.Closs(i) = out.ageing.Closs(end); %Charge lost due to side reactions
        ageing.V{i} = out.V; %Terminal voltage
        ageing.soc{i} = out.soc; %State-of-charge
        ageing.t{i} = out.t; %time vector
        
        %Display information in command window
        fprintf('Cycle %d, SoH: %f\n',i,SoH(i+1)*100)
        
        %Update variables for next simulation. The states at the end of the
        %previous simulation are stored in out.states_end, which can be
        %used as the initial condition for the next simulation.
        states_init = out.states_end;
        i = i+1; 
end
%% Plot
n_cycles = i-1;
figure()
subplot(2,1,1)
plot(0:n_cycles,SoH,'LineWidth',2)
grid on
xlabel('Cycle number [-]')
ylabel('State-of-Health [-]')

subplot(2,1,2)
for k = 1:n_cycles
    Rf(k) = mean(ageing.Rf{k}(1:out.param.grid.nn));
end
plot(1:n_cycles,Rf,'LineWidth',2)
grid on
xlabel('Cycle number [-]')
ylabel('Film resistance [Ohm m^2]')

figure()
plot(ageing.t{2}/3600,ageing.V{2},'k','LineWidth',2)
hold on
plot(ageing.t{118}/3600,ageing.V{118},'r','LineWidth',2)
plot(ageing.t{268}/3600,ageing.V{268},'b','LineWidth',2)
plot(ageing.t{470}/3600,ageing.V{470},'g','LineWidth',2)
plot(ageing.t{end}/3600,ageing.V{end},'m','LineWidth',2)
grid on
xlabel('Time [h]')
ylabel('Voltage [V]')
legend('SoH = 1', 'SoH = 0.95', 'SoH = 0.9', 'SoH = 0.85', 'SoH = 0.8')

function [i_app_next,mem,end_simulation,dt] = CC_CV(k,t,i_app,V,soc,cs,ce,phis,phie,mem,p)
V_max =p.OCV_max;   %Voltage at which the battery switches to CV mode
V_min =p.OCV_min;   %Voltage after which the battery stops discharging
i_min = 1;  %Current level to terminate the CV mode with
Kp1 = 50;   %Proportional gain of the PI controller for the CV stage
Ki1 = 10;   %Integral gain for the PI controller

end_simulation=0; %flag to indicate whether the simulation should end after the iteration. Has to be set at every iteration

% Initialize the various parameters. Note how the mem struct is used to
% store the states of the PI controller 
if k==1
    mem.status = 0;
    mem.e = 0;
    mem.e2 = 0;
    mem.end_charge = 0;
    mem.discharge_end = 1e6;
    mem.charge_end = 1e6; 
end

% Condition to terminate the initial CC discharge stage and start the CC
% charge stage
if V(k)<=V_min && mem.status==0
    mem.status = 1;  
    mem.discharge_end = k;
end

% Condition to terminate the CC charge stage and start the CV stage
if V(k)>=V_max &&mem.status==1
    mem.status = 2; 
end

%Condition to terminate the simulation
if i_app(k)<i_min && mem.status==2
    mem.charge_end = k; 
    end_simulation=1;
end

%Switch between the different stages
switch mem.status
    case 0 %CC discharge
        i_app_next = -p.discharge_current; %value of i_app for the next iteration
        dt = p.dt_discharge; %assign the chosen step size
    case 1 %CC charge
        i_app_next = p.charge_current; %value of i_app for the next iteration
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