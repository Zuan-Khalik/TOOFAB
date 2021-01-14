%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to simulate the discharge curves as shown in Fig. 7
% in [1] for the high-energy (HE) [2] cell parameters
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
addpath('Functions')
clear all; close all 
load i_app_validation.mat

soc_init = 1; 
Cap = 29.5;  

Crates = [0.1 0.5 1]; 

%% Simulate for C-rates
for k = 1:length(Crates)
    Crate = Crates(k); 
    input_current = -Cap*Crate;

    %% CDFN
    p = parameters_LS([22 10 23 3 7]); 
    p.set_simp = [1 1 1 1 0 0];
    p.dt = 1/Crate;
    out_CDFN= DFN(input_current,1e6,soc_init,p);
    V{1,k} = out_CDFN.V; 
    Q{1,k} = out_CDFN.soc*p.Cbat/3600;
    Q{1,k} = abs(Q{1,k}-Q{1,k}(1));
    ce{1,k} = out_CDFN.ce; 
    x{1,k} = out_CDFN.x; 
    jn{1,k} = out_CDFN.eta;
    t{1,k} = out_CDFN.t; 
    
    %% SDFN-HIFI
    p = parameters_LS([12 8 15 3 3]); 
    p.set_simp = [2 2 2 1 0 0];
    p.dt = 1/Crate;
    out_SDFN= DFN(input_current,1e6,soc_init,p);
    V{2,k} = out_SDFN.V; 
    Q{2,k} = out_SDFN.soc*p.Cbat/3600;
    Q{2,k} = abs(Q{2,k}-Q{2,k}(1));
    ce{2,k} = out_SDFN.ce; 
    x{2,k} = out_SDFN.x;
    jn{2,k} = out_SDFN.eta;
    t{2,k} = out_SDFN.t; 
    
    %% SPM-HIFI
    p = parameters_LS([12 8 15 3 3]); 
    p.set_simp = [2 2 2 1 0 0];
    p.dt = 1/Crate;
    out_SPM= SPM(input_current,1e6,soc_init,p);
    V{3,k} = out_SPM.V; 
    Q{3,k} = out_SPM.Q'; 
    ce{3,k} = out_SPM.ce'; 
    x{3,k} = out_SPM.x; 
    jn{3,k} = out_SPM.eta';
    t{3,k} = out_SPM.t; 
end

%% plot
plot(Q{1,1},V{1,1},'k','LineWidth',2)
hold on
plot(Q{2,1},V{2,1},'k--','LineWidth',2)
plot(Q{3,1},V{3,1},'k:','LineWidth',2)
plot(Q{1,2},V{1,2},'r','LineWidth',2)
plot(Q{2,2},V{2,2},'r--','LineWidth',2)
plot(Q{3,2},V{3,2},'r:','LineWidth',2)
plot(Q{1,3},V{1,3},'b','LineWidth',2)
plot(Q{2,3},V{2,3},'b--','LineWidth',2)
plot(Q{3,3},V{3,3},'b:','LineWidth',2)

figure()
colors{1} = {'k','k--','k:'}; 
colors{2} = {'r','r--','r:'};
colors{3} = {'b','b--','b:'};
fontsize = 16;
for k = 1:size(V,1)
    for i = 1:size(V,2)
        plot(x{k,i},ce{k,i}(:,round(0.8*length(Q{k,i})))/1200,colors{i}{k},'LineWidth',2)
        hold on
    end
end
grid on
ylabel('$c_e/c_{e,0} \ \mathrm{[-]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
save('Data/discharge_data_HE','V','Q','ce','x','jn','t')