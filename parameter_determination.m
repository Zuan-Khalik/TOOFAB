%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% V0.5.3
%
% Estimate DFN parameters
% 
% Inputs: 
%    - Equil is a struct that contains the information on the equilbrium 
%      potential of the cell (or the EMF). When the EMF is used, there should be 
%      2 fields in equil: EMF and Cbat. EMF should be function of SOC and Cbat
%      should be given in C. 
%    - data is a cell that can contain multiple structs that in itself contain
%      the measurement data. In each of the structs, there should be at least 3
%      fields: t, I, V. t is the time vector, I contains the applied current
%      measurements and V the terminal voltage measurements. For simplicity, make
%      sure that you interpolate the data such that the sample time is 1 s (see
%      the example data). 
%    - in options you can specify additional options for the parameter
%      estimation routine. In principle, you shouldn't need this field other than
%      specifying the ambient temperature. 
% 
% Outputs:
%    - The estimated parameters are given by p. 
%
% This file is a part of the TOOlbox for FAst Battery simulation (TOOFAB)
% Github: https://github.com/Zuan-Khalik/TOOFAB
%
% Author: Zuan Khalik (z.khalik@tue.nl)
%
% TOOFAB is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, ph, results] = parameter_determination(equil_data,est_data,options_in) 
options.U_neg = @(theta_n) 0.7222 + 0.1387*theta_n + 0.029*theta_n.^0.5 - 0.0172./theta_n + 0.0019./...
                        theta_n.^1.5 + 0.2808*exp(0.9-15*theta_n)-0.7984*exp(0.4465*theta_n - 0.4108);
if isfield(equil_data,'EMF') && not(isfield(equil_data,'EMF') && isfield(equil_data,'U_neg') && isfield(equil_data,'U_pos'))
    options.equil_mode = 1;
else
    options.equil_mode = 0;
end
options.sensitivity_analysis=1;
options.parameter_estimation = 1;
options.sens_perturb = 0.1; 
options.sens_perturb_mode = 1; %-1 for only lower pertubation, 1 for only higher pertubation, and 0 for both
options.plot_sensitivity = 1;
options.par_names_plot = {'$s_{0,n}$', '$s_{0,p}$','$s_{100,n}$', '$s_{100,p}$','$\hat{D}_{s,n}$', '$\hat{D}_{s,p}$', '$\hat{D}_e$', '$\hat{p}_n$', '$\hat{p}_{sep}$', '$\hat{p}_p$', '$t_+^0$',...
'$dlnfdce$','$\hat{\sigma}_n$', '$\hat{\sigma}_p$','$\hat{\kappa}$', '$\hat{R}_{f,n}$', '$\hat{R}_{cc}$',...
'$\alpha$', '$\hat{k}_{0,n}$', '$\hat{k}_{0,p}$', '$\hat{\varepsilon}_{e,n}$', '$\hat{\varepsilon}_{e,sep}$', '$\hat{\varepsilon}_{e,p}$',...
'$h_c$', '$C_p$', '$dU/dT$'}; 
options.beta_0 = 0.5;
options.rand_init =0; 
options.n_estimation = 0;
options.use_input_par_as_init = 1;
options.V_objective = 1;
options.T_objective = 0;
options.range_scaling = 0;
options.opt_tolX = 1e-4; 

options.dfn.thermal_dynamics = 0;
options.dfn.set_simp = [2 2 2 2 1 0]; 
options.dfn.T_amb = 273.15+26; 
options.dfn.verbose = 0;
options.dfn.set_grid = [10 10 10 10 10]; 

if nargin>2
    if isfield(options_in,'dfn')
        f = fieldnames(options_in.dfn);
        for i = 1:length(f)
            options.dfn.(f{i}) = options_in.dfn.(f{i});
        end
        options_in = rmfield(options_in,'dfn'); 
    end
end

if nargin>2
    f = fieldnames(options_in);
    for i = 1:length(f)
        options.(f{i}) = options_in.(f{i});
    end
end

if options.T_objective  
    options.dfn.thermal_dynamics = 1;
end

if isfield(options,'input_par')
    par_or = options.input_par; 
    out = DFN(-1,10,0.2,par_or);
    par_com = out.param; 
    options.input_par = p2phat(par_com); 
end

model = @(equil,data,par) dfn_model(data,equil,par,options); 

if options.equil_mode
par_names = {'s0_neg','s0_pos','s100_neg','s100_pos','Ds_neg', 'Ds_pos', 'De', 'p_neg', 'p_sep', 'p_pos', 't_plus',...
            'dlnfdce','sigma_neg', 'sigma_pos','kappa', 'Rf_neg', 'R_cc',...
            'alpha', 'k0_neg', 'k0_pos', 'epse_neg', 'epse_sep', 'epse_pos',...
            'hc', 'Cp', 'V_cell'}; 
else
par_names = {'Ds_neg', 'Ds_pos', 'De', 'p_neg', 'p_sep', 'p_pos', 't_plus',...
            'dlnfdce','sigma_neg', 'sigma_pos','kappa', 'Rf_neg', 'R_cc',...
            'alpha', 'k0_neg', 'k0_pos', 'epse_neg', 'epse_sep', 'epse_pos',...
            'hc', 'Cp', 'V_cell'}; 
end
ranges = determine_parameter_ranges(equil_data.Cbat);

if options.V_objective && not(options.T_objective)
    y_exp = @(exp_data) exp_data.V;
elseif options.T_objective && not(options.V_objective)
    y_exp = @(exp_data) exp_data.T;
elseif options.T_objective && options.V_objective
    y_exp = @(exp_data) [exp_data.V exp_data.T];
else 
    y_exp = @(exp_data) exp_data.V;
end
[ph,results] = estimate_parameters(model,y_exp,equil_data,par_names,ranges,est_data,options);

p = phat2p(ph,equil_data.Cbat); 
p.Cbat = equil_data.Cbat; 
if options.equil_mode
    p.EMF = equil_data.EMF;
    p.s0_neg = ph.s0_neg; 
    p.s0_pos = ph.s0_pos;
    p.s100_neg = ph.s100_neg;
    p.s100_pos = ph.s100_pos; 
else
    p.U_pos = equil_data.U_pos;
    p.U_neg = equil_data.U_neg;
    p.s0_neg = equil_data.s0_neg; 
    p.s0_pos = equil_data.s0_pos;
    p.s100_neg = equil_data.s100_neg;
    p.s100_pos = equil_data.s100_pos; 
end

f = fieldnames(options.dfn);
for i = 1:length(f)
        p.(f{i}) = options.dfn.(f{i});
end
end

function [par,results] = estimate_parameters(model,y_exp,equil_data,par_names,ranges,est_data,options)
N_par = length(par_names); %number of defined parameters of the model

%Determine the method that defines how the parameters are determined from
%beta, see Eq. (20) in [1]
for k = 1:N_par
    if ranges.([par_names{k},'_max'])/ranges.([par_names{k},'_min']) >= 10 && not(isinf(ranges.([par_names{k},'_max'])/ranges.([par_names{k},'_min'])))
        methods(k) = 1; 
    else
        methods(k) = 0;
    end
end

par = beta2par(0.5*ones(1,N_par),par_names,ranges,methods);
% Sensitivity analysis
if options.sensitivity_analysis
    dispstat(sprintf('Performing the sensitivity analysis...'),'keepthis','timestamp')
    results.sensitivity = sensitivity_analysis(model,equil_data,par_names,ranges,methods,est_data{1},options); 
end


% Parameter estimation
if options.parameter_estimation
    dispstat(sprintf('Performing the parameter estimation procedure...'),'keepthis','timestamp')
    if options.sensitivity_analysis
        if not(options.n_estimation)
        prompt = 'Please enter the number of estimation paramaters: \n'; 
        n_estimation = input(prompt); 
        else
        n_estimation = options.n_estimation; 
        end
        est_par = {results.sensitivity.ranked_parnames{1:n_estimation}}; 
        est_methods = results.sensitivity.ranked_methods; 
    else
        est_par = options.est_par; 
        n_estimation = length(est_par);
        est_methods = options.est_methods; 
    end
    beta_0 = options.beta_0*ones(1,n_estimation); %initial guess for beta
    
    %option to use a random initial guess for beta
    if options.rand_init 
        for k = 1:n_estimation
            beta_0(k) = (1+2*options.range_scaling)*rand(1)-options.range_scaling;
        end
    end

    if isfield(options,'input_par')
        f = fieldnames(options.input_par);
        for i = 1:length(f)
            par.(f{i}) = options.input_par.(f{i});
        end
        if options.use_input_par_as_init 
            for k = 1:n_estimation
                beta_0(k) = parvar_inv(par.(est_par{k}),ranges.([est_par{k},'_min']),ranges.([est_par{k},'_max']),est_methods(k)); 
                beta_0(k) = min(beta_0(k),1); 
                beta_0(k) = max(beta_0(k),0); 
            end
        end
    end
    
    eval_fun = @(beta) evaluate_objective(beta,model,par,y_exp,est_data,equil_data,est_par,est_methods,ranges,options); 
    
    options_lsqnonlin = optimoptions('lsqnonlin','Display','iter','TolFun',1e-20,'TolX',options.opt_tolX,'UseParallel',true,'MaxFunctionEvaluations',1e6,'MaxIterations',1e6);
    
    est_time_tic = tic(); 
    [beta_est,fval,~,~,lsqnonlin_output]=lsqnonlin(eval_fun,beta_0,zeros(1,n_estimation),ones(1,n_estimation),options_lsqnonlin); 
    est_time = toc(est_time_tic); 
    
    par_est = beta2par(beta_est,est_par,ranges,est_methods);
    f = fieldnames(par_est);
    for i = 1:length(f)
        par.(f{i}) = par_est.(f{i});
    end
    
    results.fval = fval;
    results.est_time = est_time; 
    results.beta_est = beta_est; 
    results.beta_0 = beta_0;
    results.lsqnonlin_output = lsqnonlin_output; 
end
end

function [F] = evaluate_objective(beta,model,par,y_exp,est_data,equil_data,est_par,est_methods,ranges,options)
par_est = beta2par(beta,est_par,ranges,est_methods);
f = fieldnames(par_est);
for i = 1:length(f)
    par.(f{i}) = par_est.(f{i});
end
 
equil = determine_equilibrium_potentials(equil_data,par,options);

p = phat2p(par,equil.Cbat); 
p.U_pos = equil.U_pos; 
p.U_neg = equil.U_neg; 
p.Cbat = equil.Cbat; 
if options.equil_mode
    p.s0_neg = par.s0_neg; 
    p.s0_pos = par.s0_pos;
    p.s100_neg = par.s100_neg;
    p.s100_pos = par.s100_pos; 
else
    p.s0_neg = equil.s0_neg; 
    p.s0_pos = equil.s0_pos;
    p.s100_neg = equil.s100_neg;
    p.s100_pos = equil.s100_pos; 
end

f = fieldnames(options.dfn);
for i = 1:length(f)
        p.(f{i}) = options.dfn.(f{i});
end

p.thermal_dynamics = 0; 

for k = 1:length(est_data)
    y_pred = model(est_data{k},equil,par); %Obtain the model output with the nominal parameters
    if length(y_pred)<length(y_exp(est_data{k})) || any(isnan(y_pred))
        y_pred = zeros(1,length(y_exp(est_data{k}))); 
    end
    
    y_exp2 = y_exp(est_data{k}); 
    Fout{k} = 10*(y_pred-y_exp2); 
end
F = Fout{1};
for k = 2:length(est_data)
    F = [F Fout{k}];
end
if any(isnan(F))
    temp = 1;
end

end

function par = beta2par(beta,par_names,ranges,method)
for i = 1:length(beta)
    if method(i)==0
        par.(par_names{i}) = (1-beta(i))*ranges.([par_names{i},'_min'])+beta(i)*ranges.([par_names{i},'_max']); 
    else
        par.(par_names{i}) = 10.^((1-beta(i))*log10(ranges.([par_names{i},'_min']))+beta(i)*log10(ranges.([par_names{i},'_max']))); 
    end
end
end

function beta = parvar_inv(theta,par_min,par_max,method)
for i = 1:length(theta)
    if method(i)==0
        beta(i) = (theta(i)-par_min(i))/(par_max(i)-par_min(i)); 
    else
        beta(i) = (log10(theta(i))-log10(par_min(i)))/(log10(par_max(i))-log10(par_min(i))); 
    end
end

end

function equil = determine_equilibrium_potentials(equil_data,par,options)
p = par;

npoints = 1000;
x = linspace(0,1,npoints); 
w = linspace(p.s0_neg, p.s100_neg,npoints);                             % negative electrode stoichiometry
z = linspace(p.s0_pos, p.s100_pos,npoints);                             % positive electrode stoichiometry

if options.equil_mode
    p.U_neg = options.U_neg; 
    p.EMF = equil_data.EMF;
    EMF = p.EMF(x); 
    U_neg = p.U_neg(w); 
    U_pos = EMF+U_neg; 
    p.EMF = griddedInterpolant(x,EMF); 
    x_2080 = x(npoints*0.2+1:npoints*0.8); 
    w_2080 = w(npoints*0.2+1:npoints*0.8); 

    temp = p.EMF(x_2080); 
    dtemp  = (temp(2:end)-temp(1:end-1))/(x(2)-x(1)); 
    EMF_mean = mean(dtemp); 

    temp = p.U_neg(w_2080); 
    dtemp  = (temp(2:end)-temp(1:end-1))/(w(2)-w(1)); 
    U_neg_mean = mean(dtemp); 

    U_pos_mean = -(EMF_mean+U_neg_mean); 

    a = U_pos_mean; 
    b = U_pos(npoints*0.2+1)-a*z(npoints*0.2+1); 

    U_pos_end = @(z) a*z+b; 

    U_pos_out1 = [U_pos_end(z(1:npoints*0.2+1)) U_pos(npoints*0.2+2:end)]; 

    U_neg_out1 = U_pos_out1-EMF; 

    U_pos_out = griddedInterpolant(flip(z),flip(U_pos_out1)); 
    U_neg_out = griddedInterpolant([0 w],[2 U_neg_out1]); 

    U_pos_full = U_pos_out(x); 
    U_neg_full = U_neg_out(x); 
    
    p.U_pos_out = griddedInterpolant(x,U_pos_full,'linear','linear'); 
    p.U_neg_out= griddedInterpolant(x,U_neg_full,'linear','linear');  
    
    p.U_pos = @(z) qinterp1(x',U_pos_full',z); 
    p.U_neg= @(w) qinterp1(x',U_neg_full',w);  
else
    U_neg = p.U_neg(w);
    U_pos = p.U_pos(z);
    EMF = U_pos-U_neg;
    p.EMF = griddedInterpolant(x,EMF);
    p.U_neg_out = p.U_neg;
    p.U_pos_out = p.U_pos;
end

equil.U_neg = p.U_neg_out; 
equil.U_pos = p.U_pos_out;
equil.Cbat = equil_data.Cbat; 

end

function results_sensitivity = sensitivity_analysis(model,equil_data,par_names,ranges,methods,sens_data,options)
%Assign the nominal parameters
N_par = length(par_names); 
beta_nom = options.beta_0*ones(N_par,1); %nominal normalized parameters
par_nom = beta2par(beta_nom,par_names,ranges,methods); %Determine the model parameter values from the normalized parameters

equil = determine_equilibrium_potentials(equil_data,par_nom,options); %determine the equilibrium potential based on the nominal parameters
y_nom = model(sens_data,equil,par_nom); %Obtain the model output with the nominal parameters

parfor k = 1:N_par
    if not(options.sens_perturb_mode==1)
    beta_perm_low = beta_nom; 
    beta_perm_low(k) = (0.5-options.sens_perturb); 
    par_perm_low = beta2par(beta_perm_low,par_names,ranges,methods);
    equil = determine_equilibrium_potentials(equil_data,par_perm_low,options); %determine the equilibrium potential based on the lower permutation of the parameters
    y_perm_low{k} = model(sens_data,equil,par_perm_low); %Obtain the model output with the lower permutation of the parameters
    dydpar_low(k,:) = -(y_perm_low{k}-y_nom)/options.sens_perturb; %Compute the sensitivity of the model output to the parameters based on the difference between the lower permutation of the parameters and the nominal parameters
    end
    
    if not(options.sens_perturb_mode==-1)
    beta_perm_high = beta_nom; 
    beta_perm_high(k) = (0.5+options.sens_perturb);
    par_perm_high = beta2par(beta_perm_high,par_names,ranges,methods);
    equil = determine_equilibrium_potentials(equil_data,par_perm_high,options); %determine the equilibrium potential based on the higher permutation of the parameters
    y_perm_high{k} = model(sens_data,equil,par_perm_high); %Obtain the model output with the higher permutation of the parameters
    dydpar_high(k,:) = (y_perm_high{k}-y_nom)/options.sens_perturb; %Compute the sensitivity of the model output to the parameters based on the difference between the higher permutation of the parameters and the nominal parameters
    end
end
switch options.sens_perturb_mode
    case -1
        S = dydpar_low; 
    case 1
        S = dydpar_high; 
    case 0
        S = 0.5*(dydpar_low+dydpar_high); %Determine the sensitivity of the model output to the parameters are the average of the two computed sensitivities 
end
[Q,R,P] = qr(S');  

Pi_init = 1:N_par; %Initial ranking of the parameters
Pi_ranked = Pi_init*P; %Compute the rank of the parameters based on the sensitivities
sensitivity_magnitude = diag(abs(R/R(1,1)));  %compute the magnitude of sensitivities

for k = 1:N_par
    ranked_parnames{k} = par_names{Pi_ranked(k)}; %par_names in ranked order
end
results_sensitivity.Pi_ranked = Pi_ranked;
results_sensitivity.Pi_init = Pi_init; 
results_sensitivity.sensitivity_magnitude = sensitivity_magnitude;
results_sensitivity.ranked_parnames = ranked_parnames; 
results_sensitivity.ranked_methods = methods(1:N_par)*P; 

if options.plot_sensitivity
    N_par_plot = sum(sensitivity_magnitude>0); 
    if isfield(options,'par_names_plot')
        par_names_plot = options.par_names_plot; 
        interpreter = 'latex'; 
        fontweight = 'bold'; 
    else
        par_names_plot = par_names; 
        interpreter = 'none'; 
        fontweight = 'normal'; 
    end

    for k = 1:N_par
        ranked_parnames_plot{k} = par_names_plot{Pi_ranked(k)};
    end
    % Plot sensitivity to output y
    figure('name','Magnitude of sensitivity')
    h = stem(sensitivity_magnitude,'k');
    set(gca,'yscal','log')
    set(gca,'ytick',10.^(-8:1:0))
    set(gca,'xtick',1:2:N_par_plot)
    xlim([0.5 N_par_plot+0.5])
    ylim([min(sensitivity_magnitude(1:N_par_plot))/10 10])
    ax = gca; 
    ax.YGrid = 'on'; 
    ax.XGrid = 'on';
    ax.YMinorGrid = 'off'; 
    h=gca; h.YAxis.TickLength = [0 0];
    set(gca,'XColor',[0 0 0])
    set(gca,'YColor',[0 0 0])

    for labelID = 1 : numel(Pi_init)
            text(Pi_init(labelID), 3*sensitivity_magnitude(labelID), ranked_parnames_plot{labelID}, 'HorizontalAlignment', 'center','Interpreter',interpreter,'FontWeight',fontweight,'Color','k');
    end
    set(gcf, 'Position',  [100, 100, 1000, 400])
    fontsize = 18;
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    xlabel('Ranking','Interpreter','latex','FontWeight','bold','Color','k','FontSize',fontsize)
    ylabel('Normalized magnitude [-]','Interpreter','latex','FontWeight','bold','Color','k','FontSize',fontsize)
    pause(1e-3)
end
end
function y = dfn_model(equil,data,ph,options)
p = phat2p(ph,equil.Cbat); 
p.U_pos = equil.U_pos; 
p.U_neg = equil.U_neg; 
p.Cbat = equil.Cbat; 
if options.equil_mode
    p.s0_neg = ph.s0_neg; 
    p.s0_pos = ph.s0_pos;
    p.s100_neg = ph.s100_neg;
    p.s100_pos = ph.s100_pos; 
else
    p.s0_neg = equil.s0_neg; 
    p.s0_pos = equil.s0_pos;
    p.s100_neg = equil.s100_neg;
    p.s100_pos = equil.s100_pos; 
end

f = fieldnames(options.dfn);
for i = 1:length(f)
        p.(f{i}) = options.dfn.(f{i});
end

p.thermal_dynamics = 0; 

if isfield(options,'Caged')
    p.ageing = 1;
    out = DFN([data.t' data.I'],data.t(end),[data.V(1), options.Caged],p); 
else
    p.ageing = 0;
    out = DFN([data.t' data.I'],data.t(end),data.V(1),p); 
end

V = interp1(out.t,out.V,data.t);  
T = interp1(out.t,out.T,data.t);  

if options.V_objective && not(options.T_objective)
   y = V; 
elseif options.T_objective && not(options.V_objective)
   y = T;
elseif options.T_objective && options.V_objective
   y = [V T];
else
    y = V;
end
    
end

function p = phat2p(ph,Cbat)
p.F = 96487; 
% Free choice parameters
rho = 20; %rho is the capacity per m^2. Taken as a typical value based on Arora et al, 2000 (rho = 23), Ecker et al, 2015 (rho = 15), Sturm et al, 2019 (rho = 26)
p.A_surf = (Cbat/3600)/rho; 
p.ce0 = 1000; 
p.delta_neg =  8.5e-5;%5e-5;
p.delta_pos = 9.2e-5;%2*3.64e-5; %6e-5;%3.64e-5; 
p.delta_sep = 2e-5; 
p.L = p.delta_sep+p.delta_neg+p.delta_pos;%2.54e-5; 
p.epss_neg = 0.53; 
p.epss_pos = 0.51;  
p.R_neg = 8.1e-6; 
p.R_pos = 5.8e-6; 

% Rest of the parameters can be determined from the grouping
p.Ds_neg = ph.Ds_neg*p.R_neg^2; 
p.Ds_pos = ph.Ds_pos*p.R_pos^2; 

p.t_plus = ph.t_plus; 
p.dlnfdce = ph.dlnfdce; 

p.De = ph.De*(1-p.t_plus)/(p.F*p.A_surf*p.ce0); 

p.epse_neg = ph.epse_neg*(1-p.t_plus)/(p.F*p.A_surf*p.delta_neg*p.ce0);
p.epse_sep = ph.epse_sep*(1-p.t_plus)/(p.F*p.A_surf*p.delta_sep*p.ce0); 
p.epse_pos = ph.epse_pos*(1-p.t_plus)/(p.F*p.A_surf*p.delta_pos*p.ce0); 

p.p_neg = log(ph.p_neg*p.delta_neg)/log(p.epse_neg); 
p.p_sep = log(ph.p_sep*p.delta_sep)/log(p.epse_sep); 
p.p_pos = log(ph.p_pos*p.delta_pos)/log(p.epse_pos);

p.sigma_neg = ph.sigma_neg*p.delta_neg/(p.epss_neg*p.A_surf); 
p.sigma_pos = ph.sigma_pos*p.delta_pos/(p.epss_pos*p.A_surf); 

p.kappa = ph.kappa/p.A_surf; 

p.Rf_neg = ph.Rf_neg*3*p.A_surf*p.delta_neg*p.epss_neg/p.R_neg; 
p.Rf_pos = 0;
p.R_cc = ph.R_cc*p.A_surf; 

p.alpha_a = ph.alpha; 
p.alpha_c = 1-p.alpha_a; 

p.k0_neg = ph.k0_neg*p.R_neg*p.F/p.ce0^p.alpha_a; 
p.k0_pos = ph.k0_pos*p.R_pos*p.F/p.ce0^p.alpha_a; 

p.V_cell = ph.V_cell; 
% p.m = 1;
p.hc = ph.hc;
p.Cp = ph.Cp;
end

function ph = p2phat(p)
p.F = 96487; 
ph.s100_neg = p.s100_neg;                                                
ph.s100_pos = p.s100_pos;                                                   
ph.s0_neg = p.s0_neg;                                                      
ph.s0_pos = p.s0_pos;    

ph.Ds_neg = p.Ds_neg/p.R_neg^2; 
ph.Ds_pos = p.Ds_pos/p.R_pos^2; 

ph.De  = p.De*p.F*p.A_surf*p.ce0/((1-p.t_plus)); 

ph.p_neg = p.epse_neg^p.p_neg/p.delta_neg; 
ph.p_sep = p.epse_sep^p.p_sep/p.delta_sep; 
ph.p_pos = p.epse_pos^p.p_pos/p.delta_pos; 

ph.t_plus = p.t_plus; 
ph.dlnfdce = p.dlnfdce; 
ph.sigma_neg = p.sigma_neg*p.epss_neg*p.A_surf/p.delta_neg; 
ph.sigma_pos = p.sigma_pos*p.epss_pos*p.A_surf/p.delta_pos; 

ph.kappa = p.kappa*p.A_surf; 

ph.Rf_neg = p.Rf_neg*p.R_neg/(3*p.A_surf*p.delta_neg*p.epss_neg); 

ph.R_cc = p.R_cc/p.A_surf;
ph.alpha = p.alpha_a ;

ph.k0_neg = p.k0_neg*p.ce0^p.alpha_a/(p.R_neg*p.F); 
ph.k0_pos = p.k0_pos*p.ce0^p.alpha_a/(p.R_pos*p.F); 

ph.epse_neg = p.epse_neg*p.F*p.A_surf*p.delta_neg*p.ce0/(1-p.t_plus); 
ph.epse_sep = p.epse_sep*p.F*p.A_surf*p.delta_sep*p.ce0/(1-p.t_plus); 
ph.epse_pos = p.epse_pos*p.F*p.A_surf*p.delta_pos*p.ce0/(1-p.t_plus); 

ph.V_cell = p.V_cell; 
ph.hc = p.hc;
ph.Cp = p.Cp;
end

function ph = determine_parameter_ranges(Cbat)
%% Redefined parameter ranges
ph.dlnfdce_min = 0; 
ph.dlnfdce_max = 5;

ph.s0_neg_min = 0.002; 
ph.s0_neg_max = 0.04; 
ph.s0_pos_min = 0.86; 
ph.s0_pos_max = 0.97; 

ph.s100_neg_min = 0.75;
ph.s100_neg_max = 0.86;
ph.s100_pos_min = 0.22; 
ph.s100_pos_max = 0.44; 

ph.Ds_neg_min = 0.00013; 
ph.Ds_neg_max = 0.0016; 
ph.Ds_pos_min = 0.0004; 
ph.Ds_pos_max = 0.63; %2; 

ph.De_min  = 4.1e-7*Cbat; %0.01*p.A_surf; 
ph.De_max  = 7.9e-6*Cbat; 

ph.p_neg_min = 50; 
ph.p_neg_max = 5700;
ph.p_sep_min = 3700; 
ph.p_sep_max = 31000;
ph.p_pos_min = 58; 
ph.p_pos_max = 4100; %7500; %4600;

ph.t_plus_min = 0.26;
ph.t_plus_max = 0.38; 

ph.sigma_neg_min = 1.2*Cbat; 
ph.sigma_neg_max =170*Cbat; 
ph.sigma_pos_min = 0.011*Cbat; %(1/100)*1000*p.A_surf; %altered
ph.sigma_pos_max = 8.3*Cbat;

ph.kappa_min = 7e-6*Cbat; 
ph.kappa_max = 3*2e-5*Cbat;

ph.Rf_neg_min = 20/Cbat; 
ph.Rf_neg_max = 330/Cbat; 

ph.R_cc_min = 32/Cbat;
ph.R_cc_max = 170/Cbat; 

ph.alpha_min = 0.48;
ph.alpha_max = 0.52 ;

ph.k0_neg_min = 5.7e-5; 
ph.k0_neg_max = 0.00078; %3*0.0005; %altered 
ph.k0_pos_min = 7.9e-5; 
ph.k0_pos_max = 5*1e-3; 

ph.epse_neg_min = 0.017*Cbat; 
ph.epse_neg_max = 0.076*Cbat; % 6*6800*p.A_surf; %altered
ph.epse_sep_min = 0.0049*Cbat; 
ph.epse_sep_max = 0.083*Cbat; 
ph.epse_pos_min = 0.01*Cbat; 
ph.epse_pos_max = 0.14*Cbat; 

ph.hc_min = 1e-3; 
ph.hc_max = 100; 
ph.Cp_min = 1;
ph.Cp_max = 5000; 
ph.V_cell_min = 1e-3;
ph.V_cell_max = 1; 
end

function dispstat(TXT,varargin)
% Prints overwritable message to the command line. If you dont want to keep
% this message, call dispstat function with option 'keepthis'. If you want to
% keep the previous message, use option 'keepprev'. First argument must be
% the message.
% IMPORTANT! In the firt call, option 'init' must be used for initialization purposes.
% Options:
%     'init'      this must be called in the begining. Otherwise, it can overwrite the previous outputs on the command line.
%     'keepthis'    the message will be persistent, wont be overwritable,
%     'keepprev'  the previous message wont be overwritten. New message will start from next line,
%     'timestamp' current time hh:mm:ss will be appended to the begining of the message.
% Example:
%   clc;
%   fprintf('12345677890\n');
%   dispstat('','init')      %Initialization. Does not print anything.
%   dispstat('Time stamp will be written over this text.'); % First output
%   dispstat('is current time.','timestamp','keepthis'); % Overwrites the previous output but this output wont be overwritten.
%   dispstat(sprintf('*********\nDeveloped by %s\n*********','Kasim')); % does not overwrites the previous output
%   dispstat('','timestamp','keepprev','keepthis'); % does not overwrites the previous output
%   dispstat('this wont be overwriten','keepthis');
%   dispstat('dummy dummy dummy');
%   dispstat('final stat');
% % Output:
%     12345677890
%     15:15:34 is current time.
%     *********
%     Developed by Kasim
%     *********
%     15:15:34 
%     this wont be overwriten
%     final stat
% Copyright (c) 2013, kasim tasdemir
% All rights reserved.


% **********
% **** Options
keepthis = 0; % option for not overwriting
keepprev = 0;
timestamp = 0; % time stamp option
init = 0; % is it initialization step?
if ~isstr(TXT)
    return
end
persistent prevCharCnt;
if isempty(prevCharCnt)
    prevCharCnt = 0;
end
if nargin == 0
    return
elseif nargin > 1
    for i = 2:nargin
        eval([varargin{i-1} '=1;']);
    end
end
if init == 1
    prevCharCnt = 0;
    return;
end
if isempty(TXT) && timestamp == 0
    return
end
if timestamp == 1
    c = clock; % [year month day hour minute seconds]
    txtTimeStamp = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
else
    txtTimeStamp = '';
end
if keepprev == 1
    prevCharCnt = 0;
end
% *************** Make safe for fprintf, replace control charachters
TXT = strrep(TXT,'%','%%');
TXT = strrep(TXT,'\','\\');
% *************** Print
TXT = [txtTimeStamp TXT '\n'];
fprintf([repmat('\b',1, prevCharCnt) TXT]);
nof_extra = length(strfind(TXT,'%%'));
nof_extra = nof_extra + length(strfind(TXT,'\\'));
nof_extra = nof_extra + length(strfind(TXT,'\n'));
prevCharCnt = length(TXT) - nof_extra; %-1 is for \n
if keepthis == 1
    prevCharCnt = 0;
end
end