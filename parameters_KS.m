%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameters of the model presented in [2] with some changes presented in
% [1]
%
% This file is a part of the TOOlbox for FAst Battery simulation (TOOFAB)
% Github: https://github.com/Zuan-Khalik/TOOFAB
%
% Github: https://github.com/Zuan-Khalik/Battery-Simulation-Toolbox
%
% Author: Zuan Khalik (z.khalik@tue.nl)
%
% TOBASIM is licensed under the BSD 3-Clause License
%
% References
% [1] Khalik et al., On trade-offs between Computational Complexity and 
% Accuracy of Electrochemistry-based Battery Models, Journal of the 
% Electrochemical Society, 2020, submitted
% [2] Smith et al., Control oriented 1d electrochemical model of lithium 
% ion battery, Energy Conversion Management, 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function p = parameters_KS(gridsize)
%-------------------------------------------------------------------------%
%--- Configuration parameters --------------------------------------------%
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%- Grid Parameters -------------------------------------------------------%
%-------------------------------------------------------------------------%
% Number of nodes (spatial discretization)
if nargin <1
    p.nn = 8;
    p.ns = 8;
    p.np = 12;
    p.nrn = 3;
    p.nrp = 5;
else
    p.nn = gridsize(1);                                                                  %In negative electrode
    p.ns = gridsize(2);                                                                   %In separator
    p.np = gridsize(3);                                                                   %In positive electrode
    p.nrn = gridsize(4);
    p.nrp = gridsize(5);
end

%-------------------------------------------------------------------------%
%- Simulation parameters -------------------------------------------------%
%-------------------------------------------------------------------------%
p.dt = 1; 
p.tol = 2e-3;                                                     %Tolerance for convergence                                                                     
p.iter_max = 1e4;                                                          %Maximum iterations for the inner loop
p.gamma = 1;                                                                %Damping coefficient for update of states
p.verbose = 1;
p.Vmin = 2.7; p.Vmax = 4.2; 
% Temperature defined to be constant for now. 
p.T_amb = 298.15;
p.T_enable = 0; 

p.fvm_method = 1; 
p.set_simp = [2 2 2 2 1 0]; 
%-------------------------------------------------------------------------%
%- Material-specific parameters ------------------------------------------%
%- Default parameters are taken from Xia et al. (2017). ------------------%
%-------------------------------------------------------------------------%
% Physical constants
p.F = 96487;                                                                %Faraday's constant [C/mol]
p.R = 8.314;                                                                %Ideal gas constant [J/mol/K]

%Lengths
p.delta_neg = 5e-5;                                                         %Neg. electrode thickness [m]
p.delta_pos = 3.64e-5;                                                      %Pos. electrode thickness [m]
p.delta_sep = 2.54e-5;                                                      %Sepator thickness [m]
p.R_neg = 1e-6;                                                               %particle radius (for both pos. and neg. electrode) [m]
p.R_pos = 1e-6;
% p.R_s = p.R_neg; 
% Porosity
p.epse_neg = 0.332;                                                         %Electrolyte volume fraction at the neg. electrode [-]
p.epse_pos = 0.33;                                                          %Electrolyte volume fraction at the pos. electrode [-]
p.epse_sep = 0.5;                                                          %Electrolyte volume fraction in the seperator [-]
p.epss_neg = 0.58;                                                            %Active material volume fraction at the neg. electrode [-]
p.epss_pos = 0.5;               
p.p_neg = 1.5; 
p.p_pos = 1.5;
p.p_sep = 1.5;    

% Diffusion coefficients
p.Ds_neg = 2e-16;                                                           %Solid-phase Li diffusion coefficient at the neg. electrode [m^2/s]
p.Ds_pos = 3.7e-16;                                                         %Solid-phase Li diffusion coefficient at the pos. electrode [m^2/s]

p.Rf0_neg = 0; 
p.Rf0_pos = 0; 
p.Rf_neg = p.Rf0_neg;
p.Rf_pos = p.Rf0_pos;
% Transport properties
p.t_plus = 0.363;                                                            %Li+ transference number [-]

% Electrode plate area
p.A_surf = 10452/100^2;                                                            %Electrode plate area [m^2]
p.R_cc = 20/100^2; 

% Transfer coefficient of surface reaction
p.alpha_a = 0.5;                                                            %Charge transfer coefficent (anodic) [-]
p.alpha_c = (1-p.alpha_a);                                                            %Charge transfer coefficient (cathodic) [-]

% Maximum concentration of Li-ion on surface of pticles
p.cs_max_neg = 16.1e3;                                                      %Maximum solid-phase concentration at the neg. electrode [mol/m^3]
p.cs_max_pos = 23.9e3;                                                      %Maximum solid-phase concentration at the pos. electrode [mol/m^3]

% Reaction rate constants
p.k0_neg =  1.38*10^(-4);                                                          %Kinetic constant in the neg. electrode [mol^(-3/2)*m^(-1/2)*s^(-1)]
p.k0_pos =   0.64*10^(-4);                                                       %Kinetic constant in the pos. electrode [mol^(-3/2)*m^(-1/2)*s^(-1)]

% Solid-phase conductivity 
p.sigma_neg = 100; 
p.sigma_pos = 10;

%Stoichiometries
p.s100_neg = 0.676;                                                   %Stoichiometry at 100% SoC at the neg. electrode [-]
p.s100_pos = 0.442;                                                     %Stoichiometry at 100% SoC at the pos. electrode [-]
p.s0_neg = 0.126;                                                      %Stoichiometry at 0% SoC at the neg. electrode [-]
p.s0_pos = 0.936;                                                      %Stoichiometry at 0% So

p.ce0 = 1200; 

p.Cap0 = (p.s100_pos-p.s0_pos)*p.epss_pos*p.delta_pos*p.A_surf*p.F*p.cs_max_pos;
p.kappa = @(c,T) 15.8e-4*c.*exp(-0.85*(1e-3*c).^1.4); 
p.dlnfdx = @(c,T) (0.601-0.24*(c/1000).^0.5+0.983.*(1-0.0052*(p.T-294))*(c/1000).^1.5)*(1-p.t_plus)^-1-1; 
p.De = @(c,T) 0.134227429226746*10e-4*10.^(-4.43-54./(p.T-229-5*(c/1000))-0.22*(c/1000));
p.Ds_pos = @(stoich,T) 55782.0157020892*10.^(-20.26+534.9*(stoich-0.5).^8+2.263*(stoich-0.5).^2); 

%Thermal parameters
p.U_T = 1.1e-3;
Vcell = 164e-3*250e-3*5e-3; 
rho_cell = 2208; 
p.m = Vcell*rho_cell; 
p.h_c = 0; 
p.C_p = 1148.5; 

p.U_neg = @(w) 8.00229 + 5.0647 * w - 12.578 * sqrt(w) - 8.6322e-4 ./ w...
        + 2.1765e-5 * ( w.^1.5 ) - 0.46016 * exp ( 15.0 * ( 0.06 - w ) )...
        - 0.55364 * exp ( -2.4326 * ( w - 0.92 ) );
p.U_pos = @(z) 85.681 * ( z.^6 ) - 357.7 * ( z.^5 ) + 613.89 * ( z.^4 ) ...
        - 555.65 * ( z.^3 ) + 281.06 * ( z.^2 ) - 76.648 * z ...
        - 0.30987 * exp ( 5.657 * ( z.^115 ) ) + 13.1983;
    
p.dU_neg=@(w)(5.0647- 12.578*(1/2)*(w).^(-1/2)+ 8.6322e-4*(w.^(-2))...
+ 2.1765e-5*(3/2) * ( w.^.5 ) + 0.46016 *15* exp ( 15.0 * ( 0.06 - w) )...
+ 0.55364* 2.4326* exp ( -2.4326 * ( w - 0.92 ) ));
p.dU_pos=@(z) (85.681*6*( z.^5 )- 357.7*5*( z.^4 )+613.89*4*( z.^3 ) ...
    - 555.65 *3* ( z.^2 ) + 281.06*2 * ( z.^1 ) - 76.648 ...
    - 0.30987 * exp ( 5.657*( z.^115 ) )*5.657*115.*(z.^114));
end
