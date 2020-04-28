%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Model Simplifications and Its Impact on Computational Complexity for an 
% Electrochemistry-Based Battery Modeling Toolbox
%
% Authors: Z. Khalik, M.C.F. Donkers, H.J. Bergveld
%
% This file is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
addpath('Functions')
addpath('Data')
clear all; close all 
load i_app_validation.mat
i_app = [1.5 0 i_app_validation(1201:end)]; 
grid_param= [10,10,10,10,10];
soc_init = 0;

Cap = 29.5; 
t_interp = [1 1251 1601:4000]; 
i_app_interp = i_app*Cap;

input_current = [t_interp' i_app_interp']; 

set_simp.full = [0 0 0 0 0 0]; 
set_simp.s1 = [0 0 0 0 1 0]; 
set_simp.s22_De = [0 2 0 0 0 0]; 
set_simp.s22_Ds = [0 0 0 2 0 0];

names = {'full','s1','s22_De','s22_Ds'}; 

p = parameters_LS(grid_param); 
for i = 1:length(names)
    p.set_simp = set_simp.(names{i}); 
    out.(names{i}) = DFN(input_current,4000,soc_init,p);
end
%%
eval_time = 3320; 

make_plots({out.full,out.s1,out.s22_De,out.s22_Ds},eval_time,1,0,0)
% 
lgd = legend({'Full model','S1','S2-$D_{e}$','S2-$D_{s,\mathrm{pos}}$'},'Interpreter','latex','FontWeight','bold'); 
lgd.NumColumns = 1;
set(gcf, 'Position',  [20, 20, 1600, 600])
set(findall(gcf,'-property','FontSize'),'FontSize',18)
lgd.FontSize = 16; 
