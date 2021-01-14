%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% V1.0
%
% Simulation of the DFN model
% 
% Inputs: 
%    - input_current: Contains information about the current profile. This 
%      field can be provided either as a scalar representing the desired 
%      applied current from time 0 to tf, an array which 
%      contains the current levels at each specified sample time, or as a
%      function which takes the output voltage, current, concentration and 
%      potentials, and the parameters as input and mainly provides the 
%      current as output. The latter form is especially useful when the
%      battery is desired to be controlled in closed-loop. Example 
%      functions for input_current are provided with the toolbox.
%    - final_time specifies the simulation time in seconds
%    - init_cond: Specifies the initial condition, which can be either an 
%      initial state-of-charge, as a value between 0 and 1, an initial 
%      voltage, or a MATLAB struct where the initial condition for 
%      a non-steady-state c_s, c_e, and T can be specified. Further 
%      details on how init_cond can be specified can be found in 
%      the documentation of the toolbox.
%    - param: Can be used to change user-configurable parameters, such as
%      all the model parameters, and simulation parameters, e.g., the
%      temporal and spatial grid discretization variables. Note that this 
%      field is optional, and a default set of parameters is already 
%      contained in the DFN function. 
% 
% Outputs:
%    - Contains all the output variables, such as the output voltage, the 
%      concentrations and the potentials.
%
% This file is a part of the TOOlbox for FAst Battery simulation (TOOFAB)
% Github: https://github.com/Zuan-Khalik/TOOFAB
%
% Author: Zuan Khalik (z.khalik@tue.nl)
%
% TOOFAB is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = DFN(input_current,final_time,init_cond,varargin)
p = default_parameters();

if nargin>3
    p = process_param(p,varargin{1}); 
end

par_or = p;
time_max = final_time; 
if isa(input_current,'function_handle')
    input_mode = 2; 
    i_app = 0; 
elseif size(input_current,2)==2
    input_mode = 1; 
    F_current = griddedInterpolant(input_current(:,1),input_current(:,2),p.current_interp,p.current_extrap); 
    i_app = F_current(0); 
elseif isscalar(input_current)
    input_mode = 0;
    i_app = input_current;
else
    error('input_current not specified correctly, please check the documentation for the proper specification of input_current')
end
% warning('on')
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:illConditionedMatrix');

%% Define variables for simulation
if nargin>3
    p = fcn_system_vectors(p,init_cond,varargin{1});
else
    p = fcn_system_vectors(p,init_cond);
end
    
par_or.EMF = p.EMF; par_or.U_pos = p.U_pos_out; par_or.U_neg = p.U_neg_out;
par_or.cs_max_neg = p.cs_max_neg; par_or.cs_max_pos = p.cs_max_pos; 
m.dummy = 0; 
m = fcn_system_matrices(p,m);

[cs_prevt, ce_prevt, phis_prev,T_prevt,Cl_prevt,Rf] = init(p,m,init_cond);
max_prealloc_size = 1e5;
phis = zeros(p.nnp,min(max_prealloc_size,round(time_max/p.dt)));
phie = zeros(p.nx,min(max_prealloc_size,round(time_max/p.dt)));
if p.set_simp(6)
    cs = zeros(p.np*2+p.nn*2,round(time_max/p.dt));
else
    zeros(p.np*p.nrp+p.nn*p.nrn,min(max_prealloc_size,round(time_max/p.dt)));
end
ce = zeros(p.nx,min(max_prealloc_size,round(time_max/p.dt)));
jn = zeros(p.nnp,min(max_prealloc_size,round(time_max/p.dt)));
eta = zeros(p.nnp,min(max_prealloc_size,round(time_max/p.dt)));
i0 = zeros(p.nnp,min(max_prealloc_size,round(time_max/p.dt)));
U = zeros(p.nnp,min(max_prealloc_size,round(time_max/p.dt)));
V = zeros(1,min(max_prealloc_size,round(time_max/p.dt)));
T = zeros(1,min(max_prealloc_size,round(time_max/p.dt)));
eta2 = zeros(p.nnp,min(max_prealloc_size,round(time_max/p.dt)));
j2 = zeros(p.nnp,min(max_prealloc_size,round(time_max/p.dt)));
if p.ageing
    se_init = [p.s100_neg; p.s100_pos; p.s0_neg; p.s0_pos]; 
    p.i02 = p.i02f(T_prevt); 
end
soc = compute_soc(cs_prevt,p); 
soc_prevt = soc; 
if soc <0 || soc>1+1e-3
    warning('Initial SOC = %d not in the range of [0,1]',soc)
end
if isa(p.k0,'function_handle')
    p.k0 = p.k0f(T_prevt); 
end
%- prev indicates the condition of the state at k-1

V(1) = phis_prev(end)-phis_prev(1);
T(1) = T_prevt; 
cs(:,1) = cs_prevt;
ce(:,1) = ce_prevt;
phis(:,1) = phis_prev;
Cl(1) = Cl_prevt; 
Rf_prevt = Rf; 
%Measure simulation time for benchmark purposes
if p.verbose>0
    dispstat('Starting simulation...','keepthis','timestamp')
end

mem.dummy = 0;
warning_set = 0;
t_vec = 0; 
end_simulation = 0; 
solution_time = 0;
num_iter = 0; 
max_iterations = 0;
dt_or = p.dt; 
%% Simulation
t= 1;
ic = 0; 
while not(end_simulation)
    %- Inner loop --------------------------------------------------------%
    tic()
%     t_vec= [t_vec t_vec(end)+p.dt]; 
    dt_prev = p.dt;
    if input_mode==2
        if nargout(input_current) == 4
            [i_app(t+1),mem,end_simulation,dt_or] = input_current(t,t_vec(t),i_app,V,soc,cs,ce,phis,phie,mem,p); 
            p.dt = dt_or;
        elseif nargout(input_current)==3
            [i_app(t+1),mem,end_simulation] = input_current(t,t_vec(t),i_app,V,soc,cs,ce,phis,phie,mem,p);
        else
            error('input_current detected as a function handle, but does not have the required inputs and outputs defined.')
        end
    end
    if max_iterations %if max iterations were exceeded in the Newton's algorithm at the previous time step
        p.dt = dt_next; %decrease the time step
    elseif not(max_iterations) && not(dt_or==p.dt)
        p.dt = dt_or; %if max iterations were not exceeded in the Newton's algorithm at the previous step (i.e., algorithm converged), choose time step as the originally chosen time step
    end
    if p.dt<1e-60 %if the time step is very small it means that the Newton's algorithm could not converge even after choosing a very small time step. 
        warning('Algorithm did not converge at time %d after lowering the sampling time. Most likely, the model is being pushed into an infeasible region (Either cs>1 or ce<0). Try changing the input or parameters. Stopping simulation',t_vec(end))
        end_simulation =1;
    end
    t_vec(t+1) = t_vec(t)+p.dt; 
    if input_mode==0
        i_app(t+1) = input_current;
    end
    if input_mode==1
        i_app(t+1) = F_current(t_vec(t+1));  
    end
    if(not(isequal(dt_prev,p.dt)))
        m = fcn_system_matrices(p,m);
    end
    
    m = fcn_system_matrices2(p,m,cs_prevt,ce_prevt); 
    
    % Solve set of AEs using Newton's method at time t
    for k=1:p.iter_max
        % Obtain function matrix and Jacobian       
        [F_phis,cs(:,t+1),ce(:,t+1),phie(:,t+1),jn(:,t+1),j2(:,t+1),i0(:,t+1),eta(:,t+1),eta2(:,t+1),U(:,t+1),p,m] = fcn_F(phis_prev,ce_prevt,cs_prevt,T_prevt,i_app(t+1),p,m);
        conv_check = norm(F_phis,2); 
        num_iter = num_iter+1; 
       % If algorithm doesn't converge, then display a warning
        if k==p.iter_max
            if not(max_iterations)
                dt_or = p.dt; 
                ic = 0; 
            end
            warning_set = 1;
            max_iterations = 1;
            ic = ic+1;
            end_simulation = 0;
            break
        end
        % If criterium for convergence is met, then break loop
        if(conv_check<p.tol ||end_simulation==1)
            max_iterations = 0;
            break
        end
        [J_phis,p] = fcn_J(cs(:,t+1),ce(:,t+1),jn(:,t+1),jn(:,t+1),i0(:,t+1),eta(:,t+1),eta2(:,t+1),T_prevt,p,m); 
        phis_prev = real(phis_prev-(J_phis\F_phis)); 
        
    end
    phis(:,t+1) = phis_prev; 
    
    V(t+1) = phis(end,t+1)+i_app(t+1)/p.A_surf*p.dx(end)*0.5/p.sigma_eff(end)...
        -phis(1,t+1)-i_app(t+1)/p.A_surf*p.dx(1)*0.5/p.sigma_eff(1)...
        +(p.R_cc/p.A_surf)*i_app(t+1);
    if mod(t,100)==0
        if p.verbose ==2
            print_loop(t_vec(end),time_max,i_app(t+1),V(t+1),k,solution_time)
        end
    end

    if p.ageing   
        Clterm = p.F*p.dt*p.dx_n*p.A_surf*sum(abs(j2(1:p.nn,t+1).*p.a_s_neg)); 
        Cl(t+1) = Cl_prevt+Clterm;
        Rf(:,t+1) = Rf_prevt-p.dt*(p.V_SEI./(p.sigma_SEI).*j2(:,t+1)); 
    end
    
    if mod(t,1000)==0 && p.ageing==1
        [Cbat_new,sbat] = fcn_Q(Cl(t+1),se_init,1,p); 
        p.s0_neg = sbat(3); p.s100_neg = sbat(1);
        p.s0_pos = sbat(4); p.s100_pos = sbat(2);
        se_init = sbat;
        Qs = soc_prevt*p.Cbat-(Cl(t+1)-Cl(t+1-1000));
        p.Cbat = Cbat_new;
        soc_prevt = Qs/p.Cbat; 
    end
    soc(t+1)= soc_prevt+p.dt*i_app(t+1)/(p.Cbat); 
    if p.set_simp(6)
        stoich = [cs(p.nn+p.np+1:p.nn+p.np+p.nn,t+1)/p.cs_max_neg; cs(2*p.nn+p.np+1:end,t+1)/p.cs_max_pos];
    else
        stoich = [cs(p.nrn:p.nrn:p.nrn*p.nn,t+1)/p.cs_max_neg; cs(p.nrn*p.nn+p.nrp:p.nrp:end,t+1)/p.cs_max_pos];
    end
    
    if p.thermal_dynamics
        T(t+1) = fcn_T(jn(:,t+1), U(:,t+1),stoich,V(t+1),i_app(t+1),T_prevt, p);
    else
        T(t+1) = p.T_amb; 
    end
    solution_time = solution_time+toc();
    if (any(stoich>=1-1e-6) || any(stoich<=1e-6)) && not(max_iterations)
        end_simulation =1;
        warning('cs exceeeding either cs_max or is lower than 0. Stopping the simulation.')
        warning_set = 1;
    end
    if any(ce(:,t+1)<=1e-6) && not(max_iterations)
        end_simulation =1;
        warning_set = 1;
        warning('ce is lower than 0. Stopping the simulation.')
    end
    if (V(t+1) <p.Vmin || V(t+1) >p.Vmax) && not(max_iterations)
        warning('Voltage exceeded specified bounds. Stopping the simulation.')
        end_simulation=1;
    end
    if end_simulation==1
        break
    end
    if t_vec(end)>=time_max
        break
    end
    
    if max_iterations
        dt_next = 0.5*p.dt; 
        if ic>3
            dt_next = 1e-20*p.dt;
        end
        phis_prev = phis(:,t);
    else
        cs_prevt = cs(:,t+1); 
        ce_prevt = ce(:,t+1); 
        T_prevt = T(t+1); 
        soc_prevt = soc(t+1); 
        if p.ageing
            Cl_prevt = Cl(t+1); 
            p.i02 = p.i02f(T_prevt); 
            Rf_prevt = Rf(:,t+1); 
            p.Rf = Rf_prevt;
        end
        if isa(p.k0,'function_handle') && p.thermal_dynamics
            p.k0 = p.k0f(T_prevt); 
        end
        t = t+1;
    end
end
out.sim_time = solution_time;
if warning_set==0
    if p.verbose>0
    dispstat(sprintf('Finished the simulation in %2.2f s \n',out.sim_time),'keepthis','timestamp');
    end
else
%     dispstat(sprintf('Finished the simulation in %2.2f s with warning %d \n',out.sim_time,warning_set),'keepthis','timestamp');
end
if warning_set ==1
    n_t = t+1;
else
    n_t = t+1; 
end
p.L = p.delta_neg+p.delta_sep+p.delta_pos; 
% Store states
out.x = [p.dx_n/2*(1:2:(p.nn*2-1)) p.delta_neg+[p.dx_s/2*(1:2:(p.ns*2-1))]...
        p.delta_neg+p.delta_sep+p.dx_p/2*(1:2:(p.np*2-1))]'/p.L; 
out.t = t_vec(2:n_t); 
out.cs = cs(:,2:n_t);
out.ce = ce(:,2:n_t); 
out.phis = [phis(1:p.nn,2:n_t); NaN(p.ns,n_t-1); phis(p.nn+1:end,2:n_t)]; 
out.phie = phie(:,2:n_t); 
if p.set_simp(6)
    out.cs_bar = [cs(p.nnp+(1:p.nn),2:n_t); NaN(p.ns,n_t-1); cs(p.nnp+(p.nn+1:p.nnp),2:n_t)];
else
    out.cs_bar = [cs(p.nrn:p.nrn:p.nrn*p.nn,2:n_t); NaN(p.ns,n_t-1); cs(p.nrn*p.nn+p.nrp:p.nrp:end,2:n_t)];
end
out.stoich = [out.cs_bar(1:p.nn,:)/p.cs_max_neg; NaN(p.ns,n_t-1); out.cs_bar(p.nns+1:end,:)/p.cs_max_pos]; 
out.jn = [jn(1:p.nn,2:n_t); NaN(p.ns,n_t-1); jn(p.nn+1:end,2:n_t)];  
out.U = [U(1:p.nn,2:n_t); NaN(p.ns,n_t-1); U(p.nn+1:end,2:n_t)];  
out.eta = [eta(1:p.nn,2:n_t); NaN(p.ns,n_t-1); eta(p.nn+1:end,2:n_t)]; 
out.V = V(2:n_t); 
out.i_app = i_app(2:n_t);
out.T = T(2:n_t); 
out.soc = soc(2:n_t);
out.param = par_or; 
out.mem = mem;

if p.ageing
    [out.ageing.Cbat_aged,sbat] = fcn_Q(Cl(n_t),se_init,1,p); 
    out.ageing.s0_neg_aged = sbat(3); 
    out.ageing.s100_neg_aged = sbat(1);
    out.ageing.s0_pos_aged = sbat(4);
    out.ageing.s100_pos_aged = sbat(2);
    out.ageing.Closs = Cl(2:n_t); 
    out.ageing.Rf = [Rf(1:p.nn,end); NaN(p.ns,1); Rf(p.nn+1:end,end)];
    out.ageing.j2 = [j2(1:p.nn,2:n_t); NaN(p.ns,n_t-1); j2(p.nn+1:end,2:n_t)]; 
    out.ageing.eta2 = [eta2(1:p.nn,2:n_t); NaN(p.ns,n_t-1); eta2(p.nn+1:end,2:n_t)]; 
    out.states_end.Closs = Cl(n_t); 
    out.states_end.Rf = out.ageing.Rf; 
end

out.states_end.phis = phis(:,n_t);
out.states_end.cs = cs(:,n_t); 
out.states_end.ce = ce(:,n_t); 
if p.thermal_dynamics
    out.states_end.T = T(:,n_t); 
end

end
%% Functions
%-------------------------------------------------------------------------%
%-- Functions for the DFN model ------------------------------------------%
%-------------------------------------------------------------------------%
function soc = compute_soc(cs,p)
ri = 1:p.nn; 

if p.set_simp(6)==1
    cs_avg_r = cs(1:p.nn); 
else
    for k = 1:p.nn
        cs_avg_r(k,:) = mean(cs(k:p.nn:p.nn*p.nrn),1); 
    end
end
cs_f = @(r)interp1([((ri-1)*p.dr_n)],cs_avg_r,r,'linear','extrap'); 

n_int = 100;
r_neg = linspace(0,p.R_neg,n_int);
dr_n = r_neg(2)-r_neg(1);
fx_kp1 = diag(r_neg(2:end).^2)*cs_f(r_neg(2:end))';
fx_k = diag(r_neg(1:end-1).^2)*cs_f(r_neg(1:end-1))'; 
cs_avg = sum(3/p.R_neg^3*(fx_kp1+fx_k)/2*dr_n);
soc = (cs_avg/p.cs_max_neg-p.s0_neg)/(p.s100_neg-p.s0_neg); 
end

function [ F_phis,cs,ce,phie,jn,j2,i0,eta,eta2,U,p,m] = fcn_F(phis,ce_prevt,cs_prevt,T,i_app,p,m)
%-------------------------------------------------------------------------%
%- Compute F_phis and J_phis ---------------------------------------------%
%-------------------------------------------------------------------------% 
jn = -m.Bphis_inv*(m.Aphis*phis+m.Cphis*i_app);
ce =  m.Gamma_ce*i_app+m.Phi_ce*phis+m.Theta_ce; 
phie = m.Gamma_phie*i_app+m.Phi_phie*phis+m.Pi_phie*log(ce); 
phie_bar = m.Aphie_bar*phie; 
if p.ageing
    eta2 = phis-phie_bar-p.U2-p.Rf.*p.F.*jn;
    j2 = -(p.i02/p.F).*exp(-2*p.alpha_c2*p.F/(p.R*T)*eta2);
    j2(p.nn+1:end) = zeros(p.np,1); 
else
    j2 = zeros(p.nn+p.np,1); 
    eta2 = zeros(p.nn+p.np,1);
end
j1 = jn-j2; 
if p.set_simp(6)==1
    cs = -m.Acs_hat_inv*(m.Bcs_hat*j1+[cs_prevt(1:p.nnp);zeros(p.nnp,1)]); 
else
    cs = -m.Acs_hat_inv*(m.Bcs_hat*j1+cs_prevt); 
end
ce_bar = m.Ace_bar*ce; 
cs_bar = m.Acs_bar*cs; 

if p.set_simp(2)==0
    m = De_matrices(ce,p.T_amb,p,m);
    m.Theta_ce = -m.Ace_hat_inv*ce_prevt; 
    m.Theta_ce_bar = m.Ace_bar*m.Theta_ce;
end

if p.set_simp(4)==0
    m = Ds_matrices(cs_bar(p.nn+1:end)/p.cs_max_pos,p.T_amb,p,m); 
    if p.set_simp(6)
        m.Theta_cs = -m.Acs_hat_inv*[cs_prevt(1:p.nnp);zeros(p.nnp,1)]; 
    else
        m.Theta_cs = -m.Acs_hat_inv*cs_prevt; 
    end
    m.Theta_cs_bar = m.Acs_bar*m.Theta_cs; 
end

if p.set_simp(1)==0
    m = kappa_matrices(ce,p.T_amb,p,m); 
    m = nu_matrices(ce,p.T_amb,p,m); 
elseif p.set_simp(3)==0
    m = nu_matrices(ce,p.T_amb,p,m);
end

i0 = p.k0.*ce_bar.^p.alpha_a.*(p.cs_bar_max-cs_bar).^p.alpha_a.*cs_bar.^p.alpha_c; 


if not(isreal(i0))
    p.warning_i0 = 1;
end

w = cs_bar(1:p.nn)/p.cs_max_neg; 
z = cs_bar(p.nn+1:end)/p.cs_max_pos; 

if not(isreal(w(1)))
    w = 0.5*ones(p.nn,1); 
    z = 0.5*ones(p.nn,1); 
end

U = [p.U_neg(w); p.U_pos(z)]; 
eta = phis-phie_bar-U-p.F*p.Rf.*jn;


if p.set_simp(5)==1 
    F_phis = p.F*j1./i0-(p.alpha_a+p.alpha_c)*p.F/(p.R*T)*eta;
else
    exp1 = exp(p.alpha_a*p.F*eta/(p.R*T)); 
    exp2 = exp(-p.alpha_c*p.F*eta/(p.R*T)); 
    F_phis = p.F*jn./i0-(exp1-exp2);
end
end

function [J_phis,p] = fcn_J(cs,ce,jn,j1,i0,eta,eta2,T,p,m)
%-------------------------------------------------------------------------%
%- Compute F_phis and J_phis ---------------------------------------------%
%-------------------------------------------------------------------------% 
%- For the sake of notation, cs_bar, ce_bar, phis, phie_bar are ----------%
%- redefined to x1,x2,x3,x4, respectively, as well as their related ------%
%- variables. ------------------------------------------------------------%

%- Exchange current density, Eq. (7) in Xia et al. (2017)

cs_bar = m.Acs_bar*cs; 
ce_bar = m.Ace_bar*ce; 
di0dcs = (-p.alpha_a./(p.cs_bar_max-cs_bar)+p.alpha_c./cs_bar).*i0; 
di0dce = p.alpha_a./ce_bar.*i0; 
djndphis = -m.Bphis_inv*m.Aphis; 

dphiedphis = m.Phi_phie_bar+m.Pi_phie_bar*diag(1./ce)*m.Phi_ce;  

if p.ageing
    deta2dphis = eye(p.nnp)-dphiedphis-diag(p.F*p.Rf)*djndphis;
    dj2dphis = diag(p.i02*2*p.alpha_c2/(p.R*T)*exp(-2*p.alpha_c2*p.F/(p.R*T)*eta2))*deta2dphis; 
    dj2dphis(p.nn+1:end,:) = zeros(p.np,p.nnp); 
    dj1dphis = djndphis-dj2dphis; 
else
    dj1dphis = djndphis; 
end

if p.ageing
    dcsdphis = -m.Acs_j1*dj1dphis;
else
    dcsdphis = m.Phi_cs_bar; 
end
dcedphis = m.Phi_ce_bar; 

di0dphis = diag(di0dcs)*dcsdphis+diag(di0dce)*dcedphis; 

w = cs_bar(1:p.nn)/p.cs_max_neg; 
z = cs_bar(p.nn+1:end)/p.cs_max_pos;

if not(isreal(w(1))) || isnan(w(1))
    w = 0.5*ones(p.nn,1); 
    z = 0.5*ones(p.nn,1); 
end

dUdcs = diag([p.dU_neg(w)/p.cs_max_neg; p.dU_pos(z)/p.cs_max_pos]); 
dUdphis = dUdcs*dcsdphis; 

detadphis = eye(p.nnp)-dphiedphis-dUdphis-diag(p.F*p.Rf)*djndphis; 

if p.set_simp(5)==1 
    J_phis = p.F*diag(1./i0)*dj1dphis-p.F*diag(j1./i0.^2)*di0dphis-(p.alpha_a+p.alpha_c)*p.F/(p.R*T)*detadphis; 
else
    exp1 = exp(p.alpha_a*p.F*eta/(p.R*T)); 
    exp2 = exp(-p.alpha_c*p.F*eta/(p.R*T)); 
    dexp1dphis = diag(p.alpha_a*p.F/(p.R*T)*exp1)*detadphis; 
    dexp2dphis = -diag(p.alpha_c*p.F/(p.R*T)*exp2)*detadphis; 
    J_phis = p.F*diag(1./i0)*djndphis-p.F*diag(jn./i0.^2)*di0dphis-(dexp1dphis-dexp2dphis); 
end
end

function [T] = fcn_T(jn, U,stoich,V,i_app,T_prevt, p)
qconvec = p.hc*(p.T_amb-T_prevt); 
DeltaT = U-T_prevt*[p.dU_dT_neg(stoich(1:p.nn)); p.dU_dT_pos(stoich(p.nn+1:end))] ; 
qchem = i_app*V-p.A_surf*p.F*ones(1,p.nnp)*(p.dx(p.elec_range).*p.a_s.*jn.*DeltaT);  %Resistive and reversible heat generation, only heat mixing due to nonuniform current distribution is accounted for
T = T_prevt+p.dt*(1/p.Cp*p.m)*(qconvec+qchem); 
end

function [ cs, ce, phis,T,Closs,Rf] = init(p,m,init_cond)
%-------------------------------------------------------------------------%
%- Initial conditions for the states cs, ce, phis, phie. -----------------%
%-------------------------------------------------------------------------%
if isstruct(init_cond)
    if not(isfield(init_cond,'cs')) || not(isfield(init_cond,'cs'))
        error('Initial condition for cs and ce are not specified in init_cond')
    else
        cs = init_cond.cs; 
        ce = init_cond.ce;
    end
    if p.thermal_dynamics
        if not(isfield(init_cond,'T'))
            error('Temperature dynamics are enabled, but initial condition for T is not specified in init_cond')
        else
            T = init_cond.T;
        end
    else
        T = p.T_amb;
    end
    phis = init_cond.phis;
    if p.ageing
        if not(isfield(init_cond,'Closs'))
            error('Ageing is enabled, but initial condition for Closs is not specified in init_cond')
        else
            Closs = init_cond.Closs;
        end
        if not(isfield(init_cond,'Rf'))
            error('Ageing is enabled, but initial condition for Rf is not specified in init_cond')
        else
            Rf = init_cond.Rf(p.elec_range);
        end
    else
        Closs = 0;
        Rf = p.Rf; 
    end
    
else
    cs0_neg = p.s0_neg*p.cs_max_neg;
    cs100_neg = p.s100_neg*p.cs_max_neg;
    cs0_pos = p.s0_pos*p.cs_max_pos;
    cs100_pos = p.s100_pos*p.cs_max_pos;

    if init_cond >2
        x = linspace(0,1,10000); 
        EMF = p.U_pos((p.s100_pos-p.s0_pos)*x+p.s0_pos)-p.U_neg((p.s100_neg-p.s0_neg)*x+p.s0_neg); 
        soc_init = interp1(EMF,x,init_cond,'linear','extrap'); 
    else
        soc_init = init_cond;
    end

    if p.set_simp(6)
        cs(1:p.nn) = (cs100_neg-cs0_neg)*soc_init+cs0_neg; 
        cs(p.nn+1:p.nnp) = (cs100_pos-cs0_pos)*soc_init+cs0_pos; 
        cs(p.nnp+1:p.nnp+p.nn) = (cs100_neg-cs0_neg)*soc_init+cs0_neg; 
        cs(p.nnp+p.nn+1:2*p.nnp) = (cs100_pos-cs0_pos)*soc_init+cs0_pos; 
        cs = cs'; 
        cs_bar= m.Acs_bar*cs;
    else

    cs(1:p.nn*p.nrn,1) = (cs100_neg-cs0_neg)*soc_init+cs0_neg;
    cs(p.nn*p.nrn+1:p.nn*p.nrn+p.np*p.nrp,1) = (cs100_pos-cs0_pos)*soc_init+cs0_pos; 
    cs_bar= [cs(p.nrn:p.nrn:p.nrn*p.nn); cs(p.nrn*p.nn+p.nrp:p.nrp:end)];
    end

    ce(:,1) = p.ce0*ones(p.nx,1);

    w = cs_bar(1:p.nn)/p.cs_max_neg; 
    z = cs_bar(p.nn+1:end)/p.cs_max_pos; 
    phis = [p.U_neg(w); p.U_pos(z)]; 

    T = p.T_amb;
    Closs = p.Closs_init;
    Rf = p.Rf; 
    
end
end

function [Q, stoich_Q] = fcn_Q(Closs,stoich_Q_prev,indicator,p)
x_prev = stoich_Q_prev; 
gamma_neg = mean(p.epss_neg)*p.delta_neg*p.cs_max_neg*p.A_surf*p.F; 
gamma_pos = mean(p.epss_pos)*p.delta_pos*p.cs_max_pos*p.A_surf*p.F;
Q0 = p.s100_neg_fresh*gamma_neg+p.s100_pos_fresh*gamma_pos; 
k = 1; 

if Closs<=1000 || indicator
    n_parts = 1;
else
    n_parts = ceil(Closs/1000);
end
Closs_orig = Closs; 
for i = 1:n_parts
    Closs = 1000*i; 
    if i==n_parts
        Closs = Closs_orig;
    end
while(1)
    F = [x_prev(1)*gamma_neg+x_prev(2)*gamma_pos-(Q0-Closs);...
        p.U_pos(x_prev(2))-p.U_neg(x_prev(1))-p.OCV_max;...
        x_prev(3)*gamma_neg+x_prev(4)*gamma_pos-(Q0-Closs);...
        p.U_pos(x_prev(4))-p.U_neg(x_prev(3))-p.OCV_min]; 
    J = [gamma_neg gamma_pos 0 0; -p.dU_neg(x_prev(1)) p.dU_pos(x_prev(2)) 0 0; ....
        0 0 gamma_neg gamma_pos;0 0 -p.dU_neg(x_prev(3)) p.dU_pos(x_prev(4))]; 
    x = x_prev-(J\F);
    if norm(x-x_prev,2)<1e-10
        break
    end
    x_prev= x;
    k = k+1; 
    if k>1000
        break
    end
end
end
stoich_Q = x; 
Q = (stoich_Q(1)-stoich_Q(3))*p.epss_neg*p.delta_neg*p.A_surf*p.cs_max_neg*p.F; 
end

function [ m ] = fcn_system_matrices( p,m )
%-------------------------------------------------------------------------%
%- This function defines the system matrices that do not change in the ---%
%- Inner loop. The matrices that do change in the inner loop are placed --%
%- in their respective functions. ----------------------------------------%
%-------------------------------------------------------------------------%
if p.set_simp(6)
    m.Acs_hat = [-eye(p.nnp) zeros(p.nnp); -eye(p.nnp) eye(p.nnp)]; 
    m.Acs_hat_inv = inv(m.Acs_hat); 
else
Acs_jn = sparse(fcn_Acs_j(p.nrn));
m.Acs_jp = sparse(fcn_Acs_j(p.nrp));
Acs_n = sparse(kron(eye(p.nn),Acs_jn)*diag((1/p.dr_n^2)*p.Ds_neg*ones(p.nn,1)'*kron(eye(p.nn), ones(1,p.nrn))));
m.Acs_hat_n = p.dt*Acs_n-speye(p.nrn*p.nn);
m.Acs_hat_inv_n = inv(m.Acs_hat_n); 

Bcs_jn = -2*(p.dr_n+p.R_neg)/(p.dr_n*p.R_neg)*[zeros(p.nrn-1,1); 1];
Bcs_jp = -2*(p.dr_p+p.R_pos)/(p.dr_p*p.R_pos)*[zeros(p.nrp-1,1); 1];
Bcs_n = sparse(kron(eye(p.nn), Bcs_jn));
Bcs_p = sparse(kron(eye(p.np), Bcs_jp));
Bcs = blkdiag(Bcs_n,Bcs_p);
m.Bcs_hat = Bcs*p.dt; 
end

Bce = (diag(p.a_s.*(1-p.t_plus)./(p.eps_e(p.elec_range)))); 
Bce = sparse([Bce(1:p.nn,:); zeros(p.ns,p.nnp); Bce(p.nn+1:end,:)]);
m.Bce_hat = p.dt*Bce; 

Aphis_neg = fcn_Aw(p.sigma_eff(1:p.nn),ones(p.nn,1),p.dx(1:p.nn),p);
Aphis_pos = fcn_Aw(p.sigma_eff(p.nn+1:end),ones(p.np,1),p.dx(p.nns+1:end),p);
m.Aphis = sparse(blkdiag(Aphis_neg, Aphis_pos)); 

Bphis = sparse(diag(-p.a_s*p.F.*p.dx(p.elec_range))); 
m.Cphis = sparse([-1/p.A_surf; zeros(p.nnp-2,1); 1/p.A_surf]); 

m.Bphie = (diag(p.a_s*p.F.*p.dx(p.elec_range))); 
m.Bphie = sparse([m.Bphie(1:p.nn,:); zeros(p.ns,p.nnp); m.Bphie(p.nn+1:end,:)]);
m.Bphie(end) = 0; 

if p.set_simp(6)
    m.Acs_bar = [zeros(p.nnp) eye(p.nnp)]; 
else
Acs_temp_n = kron(eye(p.nn),[zeros(1,p.nrn-1) 1]); 
Acs_temp_p = kron(eye(p.np),[zeros(1,p.nrp-1) 1]);
m.Acs_bar = sparse(blkdiag(Acs_temp_n,Acs_temp_p)); 
end
Ace_temp = blkdiag(eye(p.nn), zeros(p.ns), eye(p.np)); 
m.Ace_bar = sparse([Ace_temp(1:p.nn,:); Ace_temp(p.nns+1:end,:)]); 
m.Aphie_bar = m.Ace_bar;

m.Bphis_inv = inv(Bphis); 

if p.set_simp(2)==2 || p.set_simp(2)==0
    m = De_matrices(p.ce0,p.T_amb,p,m); 
end

if p.set_simp(4)==2 || p.set_simp(4)== 0
    m = Ds_matrices(ones(p.np,1),p.T_amb,p,m); 
end

if p.set_simp(1)==2 || p.set_simp(1) == 0 
    m = kappa_matrices(p.ce0,p.T_amb,p,m); 
    m = nu_matrices(p.ce0,p.T_amb,p,m); 
end
end

function [ m ] = fcn_system_matrices2( p, m,cs_prevt, ce_prevt)
%-------------------------------------------------------------------------%
if p.set_simp(2)==1
    m = De_matrices(ce_prevt,p.T_amb,p,m); 
end

if p.set_simp(4)==1
    cs_bar = m.Acs_bar*cs_prevt; 
    m = Ds_matrices(cs_bar(p.nn+1:end)/p.cs_max_pos,p.T_amb,p,m); 
end

if p.set_simp(1)==1
    m = kappa_matrices(ce_prevt,p.T_amb,p,m); 
    m = nu_matrices(ce_prevt,p.T_amb,p,m); 
elseif p.set_simp(3)==1 && not(p.set_simp(1)==0)
    m = nu_matrices(ce_prevt,p.T_amb,p,m); 
end

if p.set_simp(4)==1 || p.set_simp(4)==2 || p.set_simp(4)==0
    if p.set_simp(6)==1
        m.Theta_cs = -m.Acs_hat_inv*[cs_prevt(1:p.nnp);zeros(p.nnp,1)]; 
    else
    m.Theta_cs = -m.Acs_hat_inv*cs_prevt; 
    end
    m.Theta_cs_bar = m.Acs_bar*m.Theta_cs; 
end
if p.set_simp(2)==1 || p.set_simp(2)==2 || p.set_simp(2)==0
    m.Theta_ce = -m.Ace_hat_inv*ce_prevt; 
    m.Theta_ce_bar = m.Ace_bar*m.Theta_ce;
end
end 

function m = De_matrices(ce,T,p,m)
    m.Ace = sparse(fcn_Aw(p.De_eff(ce,T),p.eps_e.*p.dx,p.dx,p));
    m.Ace_hat = (p.dt*m.Ace-speye(p.nx)); 
    m.Ace_hat_inv = inv(m.Ace_hat);
    prem2 = m.Ace_hat_inv*(m.Bce_hat*m.Bphis_inv); 
    m.Gamma_ce = prem2*m.Cphis;
    m.Phi_ce = prem2*m.Aphis;
    m.Gamma_ce_bar = m.Ace_bar*m.Gamma_ce;
    m.Phi_ce_bar = m.Ace_bar*m.Phi_ce; 
end

function m = Ds_matrices(stoich,T,p,m)
    if p.set_simp(6)
        p.Ds = [p.Ds_neg*ones(p.nn,1); ones(p.np,1).*p.Ds_pos(stoich)]; 
        m.Bcs_hat = [-diag(3*p.dt./p.Rs); diag(p.Rs./(5*p.Ds))]; 
    else
        Acs_p = kron(speye(p.np),m.Acs_jp)*diag(((1/p.dr_p^2)*(ones(p.np,1).*p.Ds_pos(stoich,T))'*kron(speye(p.np), ones(1,p.nrp))));
        Acs_hat_p = p.dt*Acs_p-speye(p.nrp*p.np);
        Acs_hat_inv_p = inv(Acs_hat_p); 
        m.Acs_hat_inv = blkdiag(m.Acs_hat_inv_n,Acs_hat_inv_p);  
    end
    prem1 = m.Acs_hat_inv*(m.Bcs_hat*m.Bphis_inv); 
    m.Gamma_cs = prem1*m.Cphis;
    m.Phi_cs = prem1*m.Aphis;
    m.Gamma_cs_bar = m.Acs_bar*m.Gamma_cs;
    m.Phi_cs_bar = m.Acs_bar*m.Phi_cs; 
    m.Acs_j1 = m.Acs_bar*m.Acs_hat_inv*m.Bcs_hat; 
end

function m = kappa_matrices(ce,T,p,m)
    Aphie = (fcn_Aw(p.kappa_eff(ce,T),ones(p.nx,1),p.dx,p));
    Aphie(end,end-1:end) = [0 1]; 
    m.Aphie_inv = inv(Aphie); 
    prem4 = m.Aphie_inv*(m.Bphie*m.Bphis_inv); 
    m.Gamma_phie = prem4*m.Cphis;
    m.Phi_phie= prem4*m.Aphis;
    m.Gamma_phie_bar = m.Aphie_bar*m.Gamma_phie;
    m.Phi_phie_bar = m.Aphie_bar*m.Phi_phie; 
end

function m = nu_matrices(ce,T,p,m)
    Dphie = fcn_Aw(p.nu(ce,T),ones(p.nx,1),p.dx,p);
    Dphie(end,end-1:end) = [0 0]; 
    m.Pi_phie= -m.Aphie_inv*Dphie; 
%     p.Pi_phie = zeros(p.nx,p.nx); 
    m.Pi_phie_bar = m.Aphie_bar*m.Pi_phie; 
    
end

function [ Acs_j ] = fcn_Acs_j(nr )
temp1 = linspace(0,nr-2,nr-1)./linspace(1,nr-1,nr-1);
temp2 = [2 linspace(2,nr-1,nr-2)./linspace(1,nr-2,nr-2)];

temp3 = diag(temp1,-1)-2*eye(nr)+diag(temp2,1);
temp3(nr,nr-1) = 2;
Acs_j = temp3;
end
function [ Aw] = fcn_Aw(alpha1, alpha0, dx,p)
switch p.fvm_method
    case 1 %Harmonic mean
        alpha1_face= alpha1(2:end).*alpha1(1:end-1).*(dx(2:end)+dx(1:end-1))./(alpha1(2:end).*dx(1:end-1)+alpha1(1:end-1).*dx(2:end)); 
    case 2 %Linear variation
        alpha1_face = (alpha1(2:end).*dx(1:end-1)+alpha1(1:end-1).*dx(2:end))./(dx(2:end)+dx(1:end-1)); 
    case 3 %Weighted Mean
        alpha1_face = (alpha1(2:end).*dx(2:end)+alpha1(1:end-1).*dx(1:end-1))./(dx(2:end)+dx(1:end-1)); 
end

gamma = alpha1_face./(dx(2:end)+dx(1:end-1)); 

Am = diag(gamma,-1) + diag([-gamma(1); -gamma(2:end)-gamma(1:end-1); -gamma(end)])+diag(gamma,1); 

Aw0 = diag(2./(alpha0)); 
Aw = Aw0*Am; 
end
function [ p ] = fcn_system_vectors( p,init_cond,auxpar)
%-------------------------------------------------------------------------%
%- This function converts the specified parameters into vectors. ---------%
%-------------------------------------------------------------------------%
p.nn = p.grid.nn; p.ns = p.grid.ns; p.np = p.grid.np; 
p.nrn = p.grid.nrn; p.nrp = p.grid.nrp;
p.nnp = p.nn+p.np;
p.nns = p.nn+p.ns;
p.nx = p.nn + p.ns + p.np; 
p.nrt = p.nrn*p.nn+p.nrp*p.np;
p.nt = p.nrt+2*p.nx+p.nnp; 

p.dx_n = p.delta_neg / p.nn;                                                
p.dx_s = p.delta_sep / p.ns;           
p.dx_p = p.delta_pos / p.np;           

p.dr_n=p.R_neg/(p.nrn-1);
p.dr_p=p.R_pos/(p.nrp-1);


if not(p.thermal_dynamics)
    if not(isa(p.k0_neg,'function_handle')) p.k0_neg = @(T) p.k0_neg; end
    if not(isa(p.k0_pos,'function_handle')) p.k0_pos = @(T) p.k0_pos; end  
    p.k0 = [p.k0_neg(p.T_amb)*ones(p.nn,1); p.k0_pos(p.T_amb)*ones(p.np,1)];
else
    if isa(p.k0_neg,'function_handle') || isa(p.k0_pos,'function_handle')
        if not(isa(p.k0_neg,'function_handle')) p.k0_neg = @(T) p.k0_neg; end
        if not(isa(p.k0_pos,'function_handle')) p.k0_pos = @(T) p.k0_pos; end  
        p.k0 = @(T) [p.k0_neg(T)*ones(p.nn,1); p.k0_pos(T)*ones(p.np,1)];
        p.k0f = p.k0; 
    else
        p.k0 = [p.k0_neg*ones(p.nn,1); p.k0_pos*ones(p.np,1)];
    end
end

p.Rs = [p.R_neg*ones(p.nn,1); p.R_pos*ones(p.np,1)]; 
if length(p.Rf_neg)<p.nn
    p.Rf = [p.Rf_neg*ones(p.nn,1); p.Rf_pos*ones(p.np,1)];
else
    p.Rf = [p.Rf_neg; p.Rf_pos*ones(p.np,1)];
end

p.eps_e = [p.epse_neg*ones(p.nn,1); p.epse_sep*ones(p.ns,1); p.epse_pos*ones(p.np,1)];
p.eps_s = [p.epss_neg*ones(p.nn,1); p.epss_pos*ones(p.np,1)];
p.dx = [p.dx_n*ones(p.nn,1); p.dx_s*ones(p.ns,1); p.dx_p*ones(p.np,1)]; 
p.brug = [p.p_neg*ones(p.nn,1); p.p_sep*ones(p.ns,1); p.p_pos*ones(p.np,1)]; 
p.alpha_c = (1-p.alpha_a);                                                            %Charge transfer coefficient (cathodic) [-]

p.sigma_eff_neg = p.sigma_neg.*p.epss_neg;                                    %Effective solid-phase conductivity at the neg. electrode [S/m]
p.sigma_eff_pos = p.sigma_pos*p.epss_pos;                                    %Effective solid-phase conductivity at the pos. electrode [S/m]   
p.sigma_eff = [p.sigma_eff_neg*ones(p.nn,1); p.sigma_eff_pos*ones(p.np,1)];
p.a_s_neg = 3*p.epss_neg/p.R_neg;                                              %Specific interfacial surface area at the neg. electrode [1/m]
p.a_s_pos = 3*p.epss_pos/p.R_pos;                                              %Specific interfacial surface area at the pos. electrode [1/m]
p.a_s = [p.a_s_neg*ones(p.nn,1); p.a_s_pos*ones(p.np,1)];

parnames = {'kappa' 'De' 'dlnfdce' 'Ds_pos'}; 
for i = 1:length(parnames)
    if not(isa(p.(parnames{i}),'function_handle')) || p.set_simp(i)==2
        if isa(p.(parnames{i}),'function_handle')
            if i==4
                p.(parnames{i})= p.(parnames{i})((p.s100_pos+p.s0_pos)/2,p.T_amb);
            else
                p.(parnames{i}) = p.(parnames{i})(p.ce0,p.T_amb); 
            end
        end
        p.(parnames{i}) = @(c,T) p.(parnames{i}); 
        p.set_simp(i) = 2; 
    end
end

if not(isa(p.i02,'function_handle')) p.i02 = @(T) p.i02; end 

if not(isa(p.dU_dT_neg,'function_handle')) p.dU_dT_neg= @(stoich) p.dU_dT_neg*ones(p.nn,1); end
if not(isa(p.dU_dT_pos,'function_handle')) p.dU_dT_pos= @(stoich) p.dU_dT_pos*ones(p.np,1); end


p.kappa_eff = @(c,T) p.kappa(c,T).*p.eps_e.^p.brug;     %[S/m] 
p.nu = @(c,T) 2*p.R*T*p.kappa_eff(c,T)*(p.t_plus-1)/p.F.*(1+p.dlnfdce(c,T)); 
p.De_eff = @(c,T) p.De(c,T).*p.eps_e.^p.brug; 

% p.elec_range specifies the range of the electrodes throughout the cell.
p.elec_range = [linspace(1,p.nn,p.nn) linspace(p.nns+1,p.nx, p.np)];

if nargin>2
    if isfield(auxpar,'cs_max_neg')
        p.cs_max_neg = auxpar.cs_max_neg;
    end
    if isfield(auxpar,'cs_max_pos')
        p.cs_max_pos = auxpar.cs_max_pos;
    end
end

p.s100_neg_fresh = p.s100_neg; p.s0_neg_fresh = p.s0_neg; 
p.s100_pos_fresh = p.s100_pos; p.s0_pos_fresh = p.s0_pos; 

p.sd_neg = p.s100_neg_fresh-p.s0_neg_fresh; 
p.sd_pos = p.s0_pos_fresh-p.s100_pos_fresh;  
p.cs_max_neg = p.Cbat/(p.epss_neg*p.delta_neg*p.A_surf*p.F*p.sd_neg); 
p.cs_max_pos = p.Cbat/(p.epss_pos*p.delta_pos*p.A_surf*p.F*p.sd_pos); 

if nargin>2 && p.ageing==0
    if isfield(auxpar,'cs_max_neg')
        p.cs_max_neg = auxpar.cs_max_neg;
    end
    if isfield(auxpar,'cs_max_pos')
        p.cs_max_pos = auxpar.cs_max_pos;
    end
end

[p] = deter_equil(p); 


p.OCV_max = p.U_pos(p.s100_pos_fresh)-p.U_neg(p.s100_neg_fresh); 
p.OCV_min = p.U_pos(p.s0_pos_fresh)-p.U_neg(p.s0_neg_fresh); 

if p.ageing
    se_init = [p.s100_neg; p.s100_pos; p.s0_neg; p.s0_pos]; 
    if isstruct(init_cond)
        Cl_prevt = init_cond.Closs;
    else
        Cl_prevt = p.Closs_init; 
    end
    [p.Cbat,sbat] = fcn_Q(Cl_prevt,se_init,0,p);
    p.s0_neg = sbat(3); p.s100_neg = sbat(1);
    p.s0_pos = sbat(4); p.s100_pos = sbat(2);
    p.i02f = p.i02; 
end




p.cs_bar_max = [p.cs_max_neg*ones(p.nn,1); p.cs_max_pos*ones(p.np,1)];
p.cs_max = [p.cs_max_neg*ones(p.nn*p.nrn,1); p.cs_max_pos*ones(p.np*p.nrp,1)];


end

function p = deter_equil(p)
npoints = 1000;
x = linspace(0,1,npoints); 
w = linspace(p.s0_neg, p.s100_neg,npoints);                             % negative electrode stoichiometry
z = linspace(p.s0_pos, p.s100_pos,npoints);                             % positive electrode stoichiometry

if isfield(p,'EMF') && not(isfield(p,'EMF') && isfield(p,'U_neg') && isfield(p,'U_pos'))
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
    U_neg_out = griddedInterpolant(w,U_neg_out1); 

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
    U_pos_full = p.U_pos(x); 
    U_neg_full = p.U_neg(x); 
    p.U_neg_out = p.U_neg;
    p.U_pos_out = p.U_pos;
end

if size(U_neg_full,1)>1
    U_neg_full = U_neg_full';
    U_pos_full = U_pos_full';
end
dU_neg = (U_neg_full(2:end)-U_neg_full(1:end-1))./(x(2:end)-x(1:end-1)); 
dU_pos = (U_pos_full(2:end)-U_pos_full(1:end-1))./(x(2:end)-x(1:end-1)); 
p.dU_pos = @(z) qinterp1(x(2:end)',dU_pos',z); 
p.dU_neg= @(w) qinterp1(x(2:end)',dU_neg',w);  
end

function p = process_param(p,param)
parnames = {'kappa' 'De' 'dlnfdce' 'Ds_pos' 'butler_volmer' 'solid_phase_diffusion'}; 
gridnames = {'nn' 'ns' 'np' 'nrn' 'nrp'}; 

f = fieldnames(param);
for i = 1:length(f)
    if not(strcmp(f{i},'grid'))
        p.(f{i}) = param.(f{i});
    end
end

if isfield(param,'simp')
    for k = 1:6
        if isfield(param.simp,parnames{k})
            p.set_simp(k) = param.simp.(parnames{k}); 
        end
        p.simp.(parnames{k}) = p.set_simp(k); 
    end
end

if isfield(param,'set_simp')
    p.set_simp = param.set_simp; 
    for k = 1:6
        p.simp.(parnames{k}) = p.set_simp(k); 
    end
end

if isfield(param,'grid')
    for k = 1:5
        if isfield(param.grid,gridnames{k})
            p.grid.(gridnames{k}) = param.grid.(gridnames{k});
            p.set_grid(k) = p.grid.(gridnames{k}); 
        end  
    end
end

if isfield(param,'set_grid')
    for k = 1:5
        p.grid.(gridnames{k}) = param.set_grid(k); 
    end
    p.set_grid = param.set_grid; 
end

end

function p = default_parameters()
%-------------------------------------------------------------------------%
%--- Simulation parameters --------------------------------------------%
%-------------------------------------------------------------------------%
p.dt = 1;  %step size
% Number of nodes (spatial discretization)
p.grid.nn = 10; %across negative electrode
p.grid.ns = 10; %across seperator
p.grid.np = 10; %across positive electrode
p.grid.nrn = 10; %in the radial direction in the negative electrode
p.grid.nrp = 10; %in the radial direction in the positive electrode
p.set_grid = [p.grid.nn p.grid.ns p.grid.np p.grid.nrn p.grid.nrp]; %compact specification of grid parameters
 
p.tol = 1e-3; %Tolerance for convergence                                                                     
p.iter_max = 1e2; %Maximum iterations for the Newton's algorithm for solving the algebraic equations
p.Vmin = 2; p.Vmax = 5; %minimum and maximum allowable voltages. If these bounds are exceeded, the simulation stops
p.verbose = 2; %verbosity of the toolbox. Allowable settings: 0. no verbosity, 1. only indicate starting and ending of simulation with computation time, 2. in addition to 1, also show current and voltage while simulating

p.current_interp = 'linear'; %Interpolation method if input_current is given as an array. Choose any of the methods specified in the documenation of the MATLAB griddedInterpolant function
p.current_extrap = 'linear'; %Extrapolation method if input_current is given as an array. Choose any of the methods specified in the documenation of the MATLAB griddedInterpolant function

p.fvm_method = 1; %specifies how the face values of the FVM discretization are computed. Possible settings: 1: harmonic mean, 2: linear variation, 3:weighted mean.
p.thermal_dynamics = 1; % toggles thermal dynamics. Choose 0 for no thermal dynamics.
p.T_amb = 298.15; %ambient temperature
p.Closs_init = 0; %initial condition for capacity loss due to side reactions
p.ageing=0; %toggles ageing. Choose 0 for no ageing

%Specification of simplifications. See Ref 1 for how the simplifications
%are defined. Note that choosing 0 is not recommended for the 
%simplifications that can be made on the parameters kappa, De, dlnfdce,
%Ds_pos, since it negatively impacts computation time, without any 
%significant difference in model accuracy
p.simp.kappa = 1; %Choice of simplification for kappa. Choices: 0: no simplification, 1: [S2-I], 2: [S2-II]
p.simp.De = 1; %Choice of simplification for De. Choices: 0: no simplification, 1: [S2-I], 2: [S2-II]
p.simp.dlnfdce = 1; %Choice of simplification for dlnfdce. Choices: 0: no simplification, 1: [S2-I], 2: [S2-II]
p.simp.Ds_pos = 1; %Choice of simplification for Ds_pos. Choices: 0: no simplification, 1: [S2-I], 2: [S2-II]
p.simp.butler_volmer = 0; %Choice of simplification for the butler-volmer equation. Choices: 0: no simplification, 1: [S1] (linearized Butler-Volmer)
p.simp.solid_phase_diffusion = 0;

p.set_simp = [1 1 1 1 0 0]; %compact specification of grid parameters
%-------------------------------------------------------------------------%
%- DFN model parameters --------------------------------------------------%
%- Most of the default parameters have been taken from Torchio et al, 2016%
%-------------------------------------------------------------------------%
% Physical constants
p.F = 96487; %Faraday's constant [C/mol]
p.R = 8.314; %Ideal gas constant [J/mol/K]
F = p.F; R = p.R; T_amb = p.T_amb; 
%Lengths
p.delta_neg = 8.8e-5; %Neg. electrode thickness [m]
p.delta_pos = 8e-5; %Pos. electrode thickness [m]
p.delta_sep = 2.5e-5; %Separator thickness [m]           
p.R_neg = 2e-6; %particle radius in the neg. electrode [m]
p.R_pos = 2e-6;%particle radius in the pos. electrode [m]

% Porosity
p.epse_neg = 0.485; %Electrolyte volume fraction at the neg. electrode [-]
p.epse_pos = 0.385; %Electrolyte volume fraction at the pos. electrode [-]
p.epse_sep = 0.724; %Electrolyte volume fraction in the seperator [-]
p.epss_neg = 0.4824; %Active material volume fraction at the neg. electrode [-]
p.epss_pos = 0.59; %Active material volume fraction at the pos. electrode [-]              
p.p_neg = 4; %Bruggeman constant at the neg. electrode [-]
p.p_pos = 4; %Bruggeman constant at the pos. electrode [-]
p.p_sep = 4; %Bruggeman constant at the seperator [-]  

p.Ds_neg = 3.9e-14;  %Solid-phase Li diffusion coefficient at the neg. electrode [m^2/s]
%Ds_pos can be specified either as a scalar or an anonymous function with inputs concentration and temperature
p.Ds_pos = @(stoich,T) 1315383.19875943*10.^(-20.26+534.9*(stoich-0.5).^8+2.263*(stoich-0.5).^2)*exp(-5000/p.R*(1/T-1/T_amb));   %Solid-phase Li diffusion coefficient at the pos. electrode [m^2/s] 

p.Rf_neg = 5.5e-3; %SEI film resistance in the negative electrode
p.Rf_pos = 0; %SEI film resistance in the positive electrode
% Transport properties
p.t_plus = 0.364; %Li+ transference number [-]

% Electrode plate area
p.A_surf = 1; %Electrode plate area [m^2]
p.R_cc = 0; %Overall current collector resistance

% Transfer coefficient of surface reaction
p.alpha_a = 0.5; %Charge transfer coefficent (anodic) [-] Note that the cathodic transfer coefficient alpha_c is automatically computed from alpha_a+alpha_c =1

% Reaction rate constants
% Reaction rates can be specified either as a scalar or an anonymous function with input temperature
p.k0_neg = @(T) 5.031e-11*F*exp(-5000/R*(1/T-1/T_amb)); %Kinetic constant in the neg. electrode [mol^(-3/2)*m^(-1/2)*s^(-1)] Can be specified either as a scalar or an anonymous function                                                
p.k0_pos = @(T) 2.334e-11*F*exp(-5000/R*(1/T-1/T_amb)); %Kinetic constant in the pos. electrode [mol^(-3/2)*m^(-1/2)*s^(-1)] Can be specified either as a scalar or an anonymous function

% Solid-phase conductivity [S/m]
p.sigma_neg = 100; %in the negative electrode
p.sigma_pos = 100; %in the positive electrode

% Stoichiometries [-]
p.s0_neg = 0.01429; %at 0% SoC in the neg. electrode of a fresh cell
p.s100_neg = 0.85510; %at 100% SoC in the neg. electrode of a fresh cell
p.s100_pos = 0.49550; %at 0% SoC in the pos. electrode of a fresh cell
p.s0_pos = 0.99174; %at 100% SoC in the pos. electrode of a fresh cell

p.ce0 = 1000; %average electrolyte concentration [mol/m^3]

cs_max_pos = 51554; %Maximum solid-phase concentration at the pos. electrode [mol/m^3] Does not need to be specified. Just to compute Cbat

p.Cbat = -(p.s100_pos-p.s0_pos)*p.epss_pos*p.delta_pos*p.A_surf*p.F*cs_max_pos; %Reversible capacity of the battery [C]

%Transport parameters kappa, dlnfdce, De can be specified either as a 
%scalar or an anonymous function with inputs concentration and temperature
p.kappa = @(c,T)(1e-4*c.*((-10.5+0.668*1e-3*c+0.494*1e-6*c.^2) +(0.074  -1.78*1e-5*c -8.86*1e-10*c.^2).*T + (-6.96*1e-5+2.8*1e-8*c).*T.^2).^2); %Electrolyte ionic conductivity [S/m]
p.dlnfdce = @(c,T) (0.601-0.24*(c/1000).^0.5+0.983.*(1-0.0052*(T-294))*(c/1000).^1.5)*(1-p.t_plus)^-1-1; 
p.De = @(c,T) 1e-4*10.^((-4.43-54./(T-229-5e-3*c)-0.22e-3*c)); 

%Equilbrium potentials
%Must be anonymous function with input the stoichiometry 
p.U_pos = @(theta_p) (-4.656+88.669*theta_p.^2 - 401.119*theta_p.^4 + 342.909*theta_p.^6 - 462.471*theta_p.^8 + 433.434*theta_p.^10)./...
               (-1+18.933*theta_p.^2-79.532*theta_p.^4+37.311*theta_p.^6-73.083*theta_p.^8+95.96*theta_p.^10);

p.U_neg = @(theta_n) 0.7222 + 0.1387*theta_n + 0.029*theta_n.^0.5 - 0.0172./theta_n + 0.0019./...
                        theta_n.^1.5 + 0.2808*exp(0.9-15*theta_n)-0.7984*exp(0.4465*theta_n - 0.4108);
%-------------------------------------------------------------------------%
%- Lumped thermal model parameters ---------------------------------------%
%- Most of the default parameters have been taken from Torchio et al, 2016%
%-------------------------------------------------------------------------%
%Entropy change. 
%Can be either specified as a scalar value or an anonymous function with 
%input stoichiometry
s0_neg = p.s0_neg; s100_neg = p.s100_neg; 
s0_pos = p.s0_pos; s100_pos = p.s100_pos; 
p.dU_dT_neg = @(x) 0.00192209526443574*((x-s0_neg)/(s100_neg-s0_neg)).^4 + -0.00470201188333028*((x-s0_neg)/(s100_neg-s0_neg)).^3 + 0.00456069824513835*((x-s0_neg)/(s100_neg-s0_neg)).^2 + -0.00209506053870547*((x-s0_neg)/(s100_neg-s0_neg)) + 0.000210893728587236; 
p.dU_dT_pos = @(x) -0.000272434021901448*exp(0.893333198729898*((x-s100_pos)/(s0_pos-s100_pos))) + 0.000623443244648843*exp(-9.43437509677502*((x-s100_pos)/(s0_pos-s100_pos))); 

p.m = 0.45; %mass of the cell [kg]
p.hc = 0.8915; %Convective heat transfer coefficient [W/K]
p.Cp = 500; % Heat capacity [J/kg/K]

%Ageing parameters
p.U2 = 0.21; %Equilibrium potential of the side reactions
p.i02 = @(T) 10*(1.34368000000001e-11*T.^2-7.29948384000006e-09*T+9.9367796744801e-07); %Exchange current density of the side reactions
p.alpha_c2 = 0.7;
p.V_SEI = 2e-6; % Molar volume of the SEI
p.sigma_SEI = 2.3e-6; %Conductivity of the SEI
p.Cbat_init = 0;
end

%Functions for printing
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
function print_loop(t1,time_max,i_app,V,k,sim_time)
dispstat(sprintf(['-------------------------------------------\nTime: %g/',num2str(time_max),...
    ' (%.0f %%) \ni_app: %g A\nVoltage: %g V\nNumber of iterations for convergence: %d \nTime elapsed: ',num2str(sim_time,'%2.2f'),'\n-------------------------------------------'],...
    t1,t1/time_max*100,i_app,V,k));
end

function Yi = qinterp1(x,Y,xi,methodflag)
% Performs fast interpolation compared to interp1
%
% qinterp1 provides a speedup over interp1 but requires an evenly spaced
% x array.  As x and y increase in length, the run-time for interp1 increases
% linearly, but the run-time for
% qinterp1 stays constant.  For small-length x, y, and xi, qinterp1 runs about
% 6x faster than interp1.
%
%
% Usage:
%   yi = qinterp1(x,Y,xi)  - Same usage as interp1
%   yi = qinterp1(x,Y,xi,flag)
%           flag = 0       - Nearest-neighbor
%           flag = 1       - Linear (default)
%
% Example:
%   x = [-5:0.01:5];   y = exp(-x.^2/2);
%   xi = [-4.23:0.3:4.56];
%   yi = qinterp1(x,y,xi,1);
%
% Usage restrictions
%    x must be monotonically and evenly increasing
%    e.g.,  x=-36:0.02:123;
%
%    Y may be up to two-dimensional
%
% Using with non-evenly spaced arrays:
%   Frequently the user will wish to make interpolations "on the fly" from
%   a fixed pair of library (i.e., x and y) vectors.  In this case, the
%   user can generate an equally-spaced set of library data by calling
%   interp1 once, and then storing this library data in a MAT-file or
%   equivalent.  Because the speed of qinterp1 is independent of the length
%   of the library vectors, the author recommends over-sampling this
%   generated set untill memory considerations start limitting program speed.
%
%   If the user wishes to use two or more spacings (i.e., a closely-spaced
%   library in the region of fine features, and a loosely-spaced library in
%   the region of coarse features), just create multiple libraries, record
%   the switching points, and send the search data to different qinterp1
%   calls depending on its value.
%
%   Example:
%       x1 = [-5:0.01:5];   x2 = [-40:1:-5 5:1:40];
%       y1 = exp(-x1.^2/3); y2 = exp(-x2.^2/3);
%       xi = [-30:0.3:30];
%       in = xi < 5 & xi > -5;
%       yi(in) = qinterp1(x1,y1,xi(in));
%       yi(~in) = qinterp1(x2,y2,xi(~in));

% Author: N. Brahms
% Copyright 2006

% Forces vectors to be columns
x = x(:); xi = xi(:);
sx = size(x); sY = size(Y);
if sx(1)~=sY(1)
    if sx(1)==sY(2)
        Y = Y';
    else
        error('x and Y must have the same number of rows');
    end
end

if nargin>=4
    method=methodflag;
else
    method = 1;    % choose nearest-lower-neighbor, linear, etc.
                   % uses integer over string for speed
end

% Gets the x spacing
ndx = 1/(x(2)-x(1)); % one over to perform divide only once
xi = xi - x(1);      % subtract minimum of x

% Fills Yi with NaNs
s = size(Y);
if length(s)>2
    error('Y may only be one- or two-dimensional');
end
Yi = NaN*ones(length(xi),s(2));

switch method
    case 0 %nearest-neighbor method
        rxi = round(xi*ndx)+1;        % indices of nearest-neighbors
        flag = rxi<1 | rxi>length(x) | isnan(xi);
                                      % finds indices out of bounds
        nflag = ~flag;                % finds indices in bounds
        Yi(nflag,:) = Y(rxi(nflag),:);
    case 1 %linear interpolation method
        fxi = floor(xi*ndx)+1;          % indices of nearest-lower-neighbors
        flag = fxi<1 | fxi>length(x)-1 | isnan(xi);
                                        % finds indices out of bounds
        nflag = ~flag;                  % finds indices in bounds
        Yi(nflag,:) = (fxi(nflag)-xi(nflag)*ndx).*Y(fxi(nflag),:)+...
            (1-fxi(nflag)+xi(nflag)*ndx).*Y(fxi(nflag)+1,:);
end
end