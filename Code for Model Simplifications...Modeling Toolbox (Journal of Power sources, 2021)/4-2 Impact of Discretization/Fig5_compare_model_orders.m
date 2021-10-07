%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to plot Fig. 5 in [1] for the high-power (HP) [2] and the 
% high-energy (HE) cell [3]. 
%
% Model Simplifications and Their Impact on Computational Complexity for an 
% Electrochemistry-Based Battery Modeling Toolbox
%
% Authors: Z. Khalik, M.C.F. Donkers, H.J. Bergveld
%
% This file is licensed under the BSD 3-Clause License
%
% References
% [1] Khalik et al., Model Simplifications and Their Impact on Computational 
% Complexity for an Electrochemistry-Based Battery Modeling Toolbox, 
% Journal of Power Sources, 2021
% [2] Smith et al., Control oriented 1d electrochemical model of lithium 
% ion battery, Energy Conversion Management, 2007
% [3] Torchio et al., A matlab framework based on a finite novolume model
% suitable for Li-ion battery design, simulation, and control, Journal of
% the Electrochemical Society, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear all;close all;

addpath('Data')
eval_time = 950; 
fontsize = 18;
load compare_model_orders_HP
 
subplot(3,2,1)
plot(V{3},'k','LineWidth',2)
hold on
plot(V{2},'r--','LineWidth',2)
plot(V{1},'b--','LineWidth',2)
grid on
xlabel('Time [s]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
ylabel('$V_t$ [V]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
title('HP')

subplot(3,2,3)
plot(x{3},ce{3}(:,eval_time)/1200,'k','LineWidth',2)
hold on
plot(x{2},ce{2}(:,eval_time)/1200,'r--','LineWidth',2)
plot(x{1},ce{1}(:,eval_time)/1200,'b--','LineWidth',2)
ylabel('$c_e/c_{e,0}$ [-]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
grid on
set(gca,'xtick',0:0.2:1)

subplot(3,2,5)
plot(x{3},jn{3}(:,eval_time),'k','LineWidth',2)
hold on
plot(x{2},jn{2}(:,eval_time),'r--','LineWidth',2)
plot(x{1},jn{1}(:,eval_time),'b--','LineWidth',2)
xlabel('$x/L$ [-]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
ylabel('$j_n \ \mathrm{[mol/m^2/s]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
grid on
set(gca,'xtick',0:0.2:1)

NRMSE.high_HP = NRMSE_fcn(V_base,V{3})*1000; 
NRMSE.medium_HP = NRMSE_fcn(V_base,V{2})*1000; 
NRMSE.low_HP = NRMSE_fcn(V_base,V{1})*1000; 

load compare_model_orders_HE

subplot(3,2,2)
plot(V{3},'k','LineWidth',2)
hold on
plot(V{2},'r--','LineWidth',2)
plot(V{1},'b--','LineWidth',2)
grid on
xlabel('Time [s]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
title('HE')

subplot(3,2,4)
plot(x{3},ce{3}(:,eval_time)/1000,'k','LineWidth',2)
hold on
plot(x{2},ce{2}(:,eval_time)/1000,'r--','LineWidth',2)
plot(x{1},ce{1}(:,eval_time)/1000,'b--','LineWidth',2)
grid on
set(gca,'xtick',0:0.2:1)

subplot(3,2,6)
plot(x{3},jn{3}(:,eval_time),'k','LineWidth',2)
hold on
plot(x{2},jn{2}(:,eval_time),'r--','LineWidth',2)
plot(x{1},jn{1}(:,eval_time),'b--','LineWidth',2)
xlabel('$x/L$ [-]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
grid on
lgd = legend({'High','Medium','Low'},'Interpreter','latex','FontWeight','bold');
set(gca,'xtick',0:0.2:1)
title(lgd,'Model order')
lgd.NumColumns = 1;

set(gcf, 'Position',  [20, 20, 850, 950])
set(findall(gcf,'-property','FontSize'),'FontSize',18)
lgd.FontSize = 16; 
NRMSE.high_HE = NRMSE_fcn(V_base,V{3})*1000; 
NRMSE.medium_HE = NRMSE_fcn(V_base,V{2})*1000; 
NRMSE.low_HE = NRMSE_fcn(V_base,V{1})*1000; 
