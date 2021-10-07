%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Model Simplifications and Their Impact on Computational Complexity for an 
% Electrochemistry-Based Battery Modeling Toolbox
%
% Authors: Z. Khalik, M.C.F. Donkers, H.J. Bergveld
%
% This file is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear all; close all
addpath('Functions')
eta = -0.05:0.0001:0.05; 
p = parameters_LS(); 

jn_nonl = (exp(0.5*p.F/(p.R*p.T_amb)*eta)-exp(-0.5*p.F/(p.R*p.T_amb)*eta)); 
jn_lin = p.F/(p.R*p.T_amb)*eta; 

plot(eta,jn_nonl,'k','LineWidth',2)
hold on
plot(eta,jn_lin,'r--','LineWidth',2)
grid on
set(gca,'xtick',-0.05:0.025:0.05)
xlabel('$\eta \ \rm{[V]}$','Interpreter','latex','FontWeight','bold','FontSize',16)
ylabel('$j_n/i_0 \ \rm{[mol/C]}$','Interpreter','latex','FontWeight','bold','FontSize',16)

lgd = legend({'Nonlinear Butler-Volmer', 'Linearized Butler-Volmer'},'Interpreter','latex','FontWeight','bold');

set(gcf, 'Position',  [20, 20, 800, 400])
set(findall(gcf,'-property','FontSize'),'FontSize',20)

lgd.FontSize = 18; 