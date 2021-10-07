%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example showing how to use the parameter_determination function to
% estimate the DFN model parameters using experimental data. 
%
% This file is a part of the TOOlbox for FAst Battery simulation (TOOFAB)
% Github: https://github.com/Zuan-Khalik/TOOFAB
%
% Author: Zuan Khalik (z.khalik@tue.nl)
%
% TOOFAB is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


clear all; close all
addpath('Functions')
load data_cell.mat



options.dfn.T_amb = 273+25;
[p, ph,results] = parameter_determination(equil,{data{1},data{2}},options); 

%Simulate for all data sets
for k = 1:length(data)
    if not(isempty(data{k}))
        out{k} = DFN([data{k}.t' data{k}.I'],data{k}.t(end),data{k}.V(1),p); 
        V{k} = interp1(out{k}.t,out{k}.V,data{k}.t); 
        rmse(k) = rms(data{k}.V-V{k})*1000;
    end
end
ph_est = ph; 
p_est = p;
results_est = results; 

%%
close all

n_data = length(data); 
i = 1;
figure()
for k = 1:length(data)
    if not(isempty(data{k}))
        subplot(ceil(n_data/2),2,i)
        plot(data{k}.t/60,data{k}.V,'k', 'LineWidth', 2)
        hold on
        plot(out{k}.t/60,out{k}.V,'r', 'LineWidth', 2)
        i = i+1; 
        grid on
        ylabel('Voltage [V]')
        xlabel('Time [min]')
    end
end

legend('Experiment','Simulation')