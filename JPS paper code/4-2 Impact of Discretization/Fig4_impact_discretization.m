%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to plot Fig. 4 in [1] for the high-power (HP) [2] and the 
% high-energy (HE) cell [3]. 
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
% [2] Smith et al., Control oriented 1d electrochemical model of lithium 
% ion battery, Energy Conversion Management, 2007
% [3] Torchio et al., A matlab framework based on a finite novolume model
% suitable for Li-ion battery design, simulation, and control, Journal of
% the Electrochemical Society, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
addpath('Data')
clear all; close all; 

fontsize = 18;
xtick_steps = 10; 
load impact_discretization_HP

subplot(3,2,1)
plot(t,i_app/7.2,'k','LineWidth',2)
ylabel('Current [C-rate]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
xlabel('Time [s]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
grid on

subplot(3,2,3)
semilogy([grid_base{1}, grid_select],[0,rmserror{1}]*1000,'k','LineWidth',2)
hold on
semilogy([grid_base{1}, grid_select],[0,rmserror{2}]*1000,'k:','LineWidth',2)
hold on
semilogy([grid_base{1}, grid_select],[0,rmserror{3}]*1000,'k--','LineWidth',2)
hold on
semilogy([grid_base{1}, grid_select],[0,rmserror{4}]*1000,'r','LineWidth',2)
hold on
semilogy([grid_base{1}, grid_select],[0,rmserror{5}]*1000,'r--','LineWidth',2)
ylabel('NRMSE [mV]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
grid on
xlim([0 38])
title('HP')
set(gca,'ytick',10.^(-7:2:2))
set(gca,'xtick',0:xtick_steps:38)
ax = gca; 
ax.YGrid = 'on'; 
ax.XGrid = 'on';
ax.YMinorGrid = 'off'; 
h=gca; h.YAxis.TickLength = [0 0];

subplot(3,2,5)
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{1}],'k','LineWidth',2)
hold on
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{2}],'k:','LineWidth',2)
hold on
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{3}],'k--','LineWidth',2)
hold on
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{4}],'r','LineWidth',2)
hold on
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{5}],'r--','LineWidth',2)
xlabel('Value of varying grid parameter','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
ylabel('Simulation time [s]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
grid on
xlim([0 38])
set(gca,'xtick',0:xtick_steps:38)


load impact_discretization_HE
subplot(3,2,4)
semilogy([grid_base{1}, grid_select],[0,rmserror{1}]*1000,'k','LineWidth',2)
hold on
semilogy([grid_base{1}, grid_select],[0,rmserror{2}]*1000,'k:','LineWidth',2)
hold on
semilogy([grid_base{1}, grid_select],[0,rmserror{3}]*1000,'k--','LineWidth',2)
hold on
semilogy([grid_base{1}, grid_select],[0,rmserror{4}]*1000,'r','LineWidth',2)
hold on
semilogy([grid_base{1}, grid_select],[0,rmserror{5}]*1000,'r--','LineWidth',2)
% legend('Vary n_n', 'Vary n_s', 'Vary n_p', 'Vary n_{r,n}', 'Vary n_{r,p}')
legend({'Vary $n_n$','Vary $n_s$','Vary $n_p$','Vary $n_{r,n}$', 'Vary $n_{r,p}$'},'Interpreter','latex','FontWeight','bold')

% ylabel('NRMSE [mV]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
grid on
xlim([0 38])
title('HE')
set(gca,'ytick',10.^(-7:2:2))
set(gca,'xtick',0:xtick_steps:38)
ax = gca; 
ax.YGrid = 'on'; 
ax.XGrid = 'on';
ax.YMinorGrid = 'off'; 
h=gca; h.YAxis.TickLength = [0 0];

subplot(3,2,6)
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{1}],'k','LineWidth',2)
hold on
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{2}],'k:','LineWidth',2)
hold on
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{3}],'k--','LineWidth',2)
hold on
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{4}],'r','LineWidth',2)
hold on
plot([grid_base{1}, grid_select],[sim_time_base,sim_time_mean{5}],'r--','LineWidth',2)
xlabel('Value of varying grid parameter','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
% ylabel('Simulation time [s]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
grid on
xlim([0 38])
set(gca,'xtick',0:xtick_steps:38)

set(gcf, 'Position',  [20, 20, 850, 950])
set(findall(gcf,'-property','FontSize'),'FontSize',18)