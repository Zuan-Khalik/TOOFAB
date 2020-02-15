%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to plot the discharge curves as shown in Fig. 9 in [1] for the
% high-power (HP) [3] and high-energy (HE) [2] parameters
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

addpath('Functions')
addpath('Data')
clear all; close all 
load discharge_data_HP.mat
colors{1} = {'k','k--','k:'}; 
colors{2} = {'r','r--','r:'};
colors{3} = {'b','b--','b:'};
fontsize = 16;
%% plot V HP
subplot(3,2,1)
for k = 1:size(V,1)
    for i = 1:size(V,2)
        plot(Q{k,i},V{k,i},colors{i}{k},'LineWidth',2)
        hold on
    end
end
grid on
ylabel('$V_t \ \mathrm{[V]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
xlabel('Discharged capacity [Ah]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
title('HP')
%% plot ce HP
subplot(3,2,3)
for k = 1:size(V,1)
    for i = 1:size(V,2)
        plot(x{k,i},ce{k,i}(round(0.8*length(Q{k,i})),:)/1200,colors{i}{k},'LineWidth',2)
        hold on
    end
end
grid on
set(gca,'xtick',0:0.2:1)
ylabel('$c_e/c_{e,0} \ \mathrm{[-]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)

%% plot eta HP
subplot(3,2,5)
for k = 1:size(V,1)
    for i = 1:size(V,2)
        plot(x{k,i},jn{k,i}(round(0.8*length(Q{k,i})),:),colors{i}{k},'LineWidth',2)
        hold on
    end
end
grid on
set(gca,'xtick',0:0.2:1)
ylabel('$j_n \ \mathrm{[mol/m^2/s]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
xlabel('$x/L$ \ [-]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
legend({'1C','10C','20C'},'Interpreter','latex','FontWeight','bold')
%% load data HE
load discharge_data_HE.mat
%% plot V HE
subplot(3,2,2)
for k = 1:size(V,1)
    for i = 1:size(V,2)
        plot(Q{k,i},V{k,i},colors{i}{k},'LineWidth',2)
        hold on
    end
end
grid on
xlabel('Discharged capacity [Ah]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
title('HE')
%% plot ce HE
subplot(3,2,4)
for k = 1:size(V,1)
    for i = 1:size(V,2)
        plot(x{k,i},ce{k,i}(round(0.8*length(Q{k,i})),:)/1200,colors{i}{k},'LineWidth',2)
        hold on
    end
end
set(gca,'xtick',0:0.2:1)
grid on
%% plot eta HE
subplot(3,2,6)
for k = 1:size(V,1)
    for i = 1:size(V,2)
        plot(x{k,i},jn{k,i}(round(0.8*length(Q{k,i})),:),colors{i}{k},'LineWidth',2)
        hold on
    end
end
grid on
set(gca,'xtick',0:0.2:1)
xlabel('$x/L$\ [-]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
legend({'0.1C','0.5C','1C'},'Interpreter','latex','FontWeight','bold')

set(gcf, 'Position',  [20, 20, 850, 950])
set(findall(gcf,'-property','FontSize'),'FontSize',20)