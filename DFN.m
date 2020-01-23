%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulation of the DFN model
% 
% Inputs: 
%    - input_current: contains information about the current profile. 
%      This field can be provided either as an array which contains the 
%      current levels at each specified time sample, or as a function 
%      which takes the output voltage, current, concentration and 
%      potentials, and the parameters as input and mainly provides the 
%      current as output. The latter form is especially useful when the 
%      battery is desired to be controlled in closed-loop.
%    - tf specifies the simulation time in seconds
%    - init_cond specifies the initial condition, which can be either an 
%      initial State-of-Charge (SoC), as a value between 0 and 1, or an 
%      initial voltage. If init_cond is chosen to be larger than 2, the 
%      toolbox assumes that init_cond specifies an initial voltage.
%    - param is a struct containing all the model parameters, and 
%      simulation parameters, such as the temporal and spatial grid 
%      discretization variables.
% 
% Outputs:
%    - out is a MATLAB struct containing all the output variables, such as 
%      the output voltage, the concentrations and the potentials.
%
% This file is a part of the TOOlbox for BAttery SIMulation (TOBASIM)
% Github: https://github.com/Zuan-Khalik/Battery-Simulation-Toolbox
%
% Author: Zuan Khalik (z.khalik@tue.nl)
%
% TOBASIM is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function out = DFN(input_current,tf,init_cond,param)
p = param; 
time_max = tf; 
if isa(input_current,'function_handle')
    input_mode = 2; 
    i_app = 0; 
else
    input_mode = 1; 
    F_current = griddedInterpolant(input_current(:,1),input_current(:,2),'previous','nearest'); 
    i_app = F_current(p.dt); 
end
% warning('on')
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:illConditionedMatrix');

%% Define variables for simulation
soc = init_cond;
p = fcn_system_vectors(p);
p = fcn_system_matrices(p);

[cs_prevt, ce_prevt, phis_prev,T_prevt,soc_prevt] = init(p,init_cond);
phis = zeros(p.nnp,time_max);
phie = zeros(p.nx,time_max);
if p.set_simp(6)
    cs = zeros(p.np*2+p.nn*2,time_max);
else
    cs = zeros(p.np*p.nrp+p.nn*p.nrn,time_max);
end
ce = zeros(p.nx,time_max);
jn = zeros(p.nnp,time_max);
eta = zeros(p.nnp,time_max);
U = zeros(p.nnp,time_max);
V = zeros(1,time_max);
T = zeros(1,time_max);
%- prev indicates the condition of the state at k-1

%Measure simulation time for benchmark purposes
if p.verbose<2
    dispstat('Starting simulation...','keepthis','timestamp')
end

warning_set = 0;
t_vec = p.dt; 
end_simulation = 0; 
solution_time = 0;
num_iter = 0; 
%% Simulation
for t=1:1e6
    t1=t*p.dt;
    %- Inner loop --------------------------------------------------------%
    tic()
    p = fcn_system_matrices2(p,cs_prevt,ce_prevt); 
    % Solve set of AEs using Newton's method at time t
    for k=1:p.iter_max
        % Obtain function matrix and Jacobian
        cs(:,t) = p.Gamma_cs*i_app(t)+p.Phi_cs*phis_prev+p.Theta_cs; 
        ce(:,t) = p.Gamma_ce*i_app(t)+p.Phi_ce*phis_prev+p.Theta_ce; 
        phie(:,t) = p.Gamma_phie*i_app(t)+p.Phi_phie*phis_prev+p.Pi_phie*log(ce(:,t)); 
        [F_phis,jn(:,t),i0,eta(:,t),U(:,t),p] = fcn_F(phis_prev,cs(:,t),ce(:,t),phie(:,t),ce_prevt,cs_prevt,T_prevt, i_app(t),p);
        conv_check = norm(F_phis,2); 
        num_iter = num_iter+1; 
       % If algorithm doesn't converge, then display a warning
        if k==p.iter_max
            warning_set = 1;
            end_simulation = 1;
        end 
        % If criterium for convergence is met, then break loop
        if(conv_check<p.tol ||end_simulation==1)
            if mod(t1,100)==0
                if p.verbose ==0
                print_loop(t1,time_max,i_app(t),k,solution_time)
                end
            end
            break
        end
        [J_phis,p] = fcn_J(cs(:,t),ce(:,t),jn(:,t),i0,eta(:,t), ce_prevt,T_prevt,p); 
        phis_prev = phis_prev-(J_phis\F_phis); 
        
    end
    phis(:,t) = phis_prev; 

    solution_time = solution_time+toc();
    V(t) = phis(end,t)+i_app(t)/p.A_surf*p.dx(end)*0.5/p.sigma_eff(end)...
        -phis(1,t)-i_app(t)/p.A_surf*p.dx(1)*0.5/p.sigma_eff(1)...
        +(p.R_cc/p.A_surf)*i_app(t);
    if p.T_enable
        T(t) = fcn_T(jn(:,t), U(:,t),V(t),i_app(t),T_prevt, p);
    else
        T(t) = p.T_amb; 
    end

    soc(t)= soc_prevt+p.dt*i_app(t)/(p.Cap0); 
    t_vec= [t_vec t_vec(end)+p.dt]; 
    
    if input_mode==2
        [i_app(t+1),end_simulation,p] = input_current(t1,V,soc,i_app,cs,ce,phis,phie,p); 
    else
    if not(t1==time_max+p.dt)
        i_app(t+1) = F_current(t1+p.dt);  
    end
    end
    dt_prev = t_vec(t+1)-t_vec(t); 
    if(not(isequal(dt_prev,p.dt)))
        p = fcn_system_matrices(p);
    end
    if V(t) <p.Vmin || V(t) >p.Vmax
        end_simulation=1;
    end
    if end_simulation==1
        break
    end
    if t_vec(end)>=time_max+p.dt
        break
    end
    cs_prevt = cs(:,t); 
    ce_prevt = ce(:,t); 
    T_prevt = T(t); 
    soc_prevt = soc(t); 
end
out.sim_time = solution_time;
out.num_iter = num_iter; 
if warning_set==0
    if p.verbose<2
    dispstat(sprintf('Finished the simulation in %2.2f s \n',out.sim_time),'keepthis','timestamp');
    end
else
    dispstat(sprintf('Finished the simulation in %2.2f s with warning %d \n',out.sim_time,warning_set),'keepthis','timestamp');
end
i_app = i_app(1:end-1);
n_t = t; 
p.L = p.delta_neg+p.delta_sep+p.delta_pos; 
% Store states
out.x = [p.dx_n/2*(1:2:(p.nn*2-1)) p.delta_neg+[p.dx_s/2*(1:2:(p.ns*2-1))]...
        p.delta_neg+p.delta_sep+p.dx_p/2*(1:2:(p.np*2-1))]/p.L; 
out.phis = [phis(1:p.nn,1:n_t); NaN(p.ns,n_t); phis(p.nn+1:end,1:n_t)]'; 
out.ce = ce(:,1:n_t)'; 
out.cs = cs(:,1:n_t)';  
out.phie = phie(:,1:n_t)'; 
if p.set_simp(6)
    out.cs_bar = [cs(p.nnp+(1:p.nn),1:n_t); NaN(p.ns,n_t); cs(p.nnp+(p.nn+1:p.nnp),1:n_t)]';
else
    out.cs_bar = [cs(p.nrn:p.nrn:p.nrn*p.nn,1:n_t); NaN(p.ns,n_t); cs(p.nrn*p.nn+p.nrp:p.nrp:end,1:n_t)]';
end
out.stoich = [out.cs_bar(:,1:p.nn)/p.cs_max_neg NaN(n_t,p.ns) out.cs_bar(:,p.nns+1:end)/p.cs_max_pos]; 
out.jn = [jn(1:p.nn,1:n_t); NaN(p.ns,n_t); jn(p.nn+1:end,1:n_t)]';  
out.U = [U(1:p.nn,1:n_t); NaN(p.ns,n_t); U(p.nn+1:end,1:n_t)]';  
out.eta = [eta(1:p.nn,1:n_t); NaN(p.ns,n_t); eta(p.nn+1:end,1:n_t)]'; 
out.V = V(1:n_t)'; 
out.i_app = i_app';
out.T = T(1:n_t)'; 
out.p = p; 
out.t = t_vec(1:end-1)'; 
out.soc = soc';
out.Q = out.soc*p.Cap0/3600;
out.Q = abs(out.Q-out.Q(1));
end
%% Functions
%-------------------------------------------------------------------------%
%-- Functions for the DFN model ------------------------------------------%
%-------------------------------------------------------------------------%
function [ F_phis,jn,i0,eta,U,p] = fcn_F(phis,cs,ce,phie, ce_prevt,cs_prevt,T, i_app,p)
%-------------------------------------------------------------------------%
%- Compute F_phis and J_phis ---------------------------------------------%
%-------------------------------------------------------------------------% 

cs = p.Acs_bar*cs; 
phie = p.Aphie_bar*phie; 

if p.set_simp(2)==0
    p = De_matrices(ce,p.T_amb,p);
    p.Theta_ce = -p.Ace_hat_inv*ce_prevt; 
    p.Theta_ce_bar = p.Ace_bar*p.Theta_ce;
end

if p.set_simp(4)==0
    cs_bar = cs; 
    p = Ds_matrices(cs_bar(p.nn+1:end)/p.cs_max_pos,p.T_amb,p); 
    if p.set_simp(6)
        p.Theta_cs = -p.Acs_hat_inv*[cs_prevt(1:p.nnp);zeros(p.nnp,1)]; 
    else
        p.Theta_cs = -p.Acs_hat_inv*cs_prevt; 
    end
    p.Theta_cs_bar = p.Acs_bar*p.Theta_cs; 
end

if p.set_simp(1)==0
    p = kappa_matrices(ce,p.T_amb,p); 
    p = nu_matrices(ce,p.T_amb,p); 
elseif p.set_simp(3)==0
    p = nu_matrices(ce,p.T_amb,p);
end
ce = p.Ace_bar*ce; 
%- Exchange current density, Eq. (7) in Xia et al. (2017)

jn = -p.Bphis_inv*(p.Aphis*phis+p.Cphis*i_app); 

i0 = p.k0.*ce.^p.alpha_a.*(p.cs_bar_max-cs).^p.alpha_a.*cs.^p.alpha_c; 

w = cs(1:p.nn)/p.cs_max_neg; 
z = cs(p.nn+1:end)/p.cs_max_pos; 
U = [p.U_neg(w); p.U_pos(z)]; 
eta = phis-phie-U-p.F*p.R_f.*jn;

if p.set_simp(5)==1 
    F_phis = p.F*jn./i0-(p.alpha_a+p.alpha_c)*p.F/(p.R*T)*eta;
else
    exp1 = exp(p.alpha_a*p.F*eta/(p.R*T)); 
    exp2 = exp(-p.alpha_c*p.F*eta/(p.R*T)); 
    F_phis = p.F*jn./i0-(exp1-exp2);
end
end

function [J_phis,p] = fcn_J(cs,ce,jn,i0,eta, ce_prevt,T,p)
%-------------------------------------------------------------------------%
%- Compute F_phis and J_phis ---------------------------------------------%
%-------------------------------------------------------------------------% 
%- For the sake of notation, cs_bar, ce_bar, phis, phie_bar are ----------%
%- redefined to x1,x2,x3,x4, respectively, as well as their related ------%
%- variables. ------------------------------------------------------------%

%- Exchange current density, Eq. (7) in Xia et al. (2017)

ce1 = ce; 
cs = p.Acs_bar*cs; 
ce = p.Ace_bar*ce; 
di0dcs = (-p.alpha_a./(p.cs_bar_max-cs)+p.alpha_c./cs).*i0; 
di0dce = p.alpha_a./ce.*i0; 
djndphis = -p.Bphis_inv*p.Aphis; 

dphiedphis = p.Phi_phie_bar+p.Pi_phie_bar*diag(1./ce1)*p.Phi_ce;  

dcsdphis = p.Phi_cs_bar; 
dcedphis = p.Phi_ce_bar; 

di0dphis = diag(di0dcs)*dcsdphis+diag(di0dce)*dcedphis; 

dUdcs = diag([p.dU_neg(cs(1:p.nn)/p.cs_max_neg)/p.cs_max_neg; p.dU_pos(cs(p.nn+1:end)/p.cs_max_pos)/p.cs_max_pos]); 
dUdphis = dUdcs*dcsdphis; 


detadphis = eye(p.nnp)-dphiedphis-dUdphis-diag(p.F*p.R_f)*djndphis; 

if p.set_simp(5)==1 
    J_phis = p.F*diag(1./i0)*djndphis-p.F*diag(jn./i0.^2)*di0dphis-(p.alpha_a+p.alpha_c)*p.F/(p.R*T)*detadphis; 
else
    exp1 = exp(p.alpha_a*p.F*eta/(p.R*T)); 
    exp2 = exp(-p.alpha_c*p.F*eta/(p.R*T)); 
    dexp1dphis = diag(p.alpha_a*p.F/(p.R*T)*exp1)*detadphis; 
    dexp2dphis = -diag(p.alpha_c*p.F/(p.R*T)*exp2)*detadphis; 
    J_phis = p.F*diag(1./i0)*djndphis-p.F*diag(jn./i0.^2)*di0dphis-(dexp1dphis-dexp2dphis); 
end
end

function [T] = fcn_T(jn, U,V,i_app,T_prevt, p)
qconvec = p.h_c*(p.T_amb-T_prevt); 
DeltaT = U-T_prevt*p.U_T; 
qrev = i_app*V; 
qchem = p.A_surf*p.F*ones(1,p.nnp)*(p.dx(p.elec_range).*p.a_s.*jn.*DeltaT); 
T = T_prevt+p.dt*(1/p.C_p*p.m)*(qconvec+qrev-qchem); 
end

function [ cs, ce, phis,T,soc] = init(p,init_cond)
%-------------------------------------------------------------------------%
%- Initial conditions for the states cs, ce, phis, phie. -----------------%
%-------------------------------------------------------------------------%
cs0_neg = p.s0_neg*p.cs_max_neg;
cs100_neg = p.s100_neg*p.cs_max_neg;
cs0_pos = p.s0_pos*p.cs_max_pos;
cs100_pos = p.s100_pos*p.cs_max_pos;

if init_cond >2
    x = linspace(0,1,10000); 
    EMF = p.U_pos((p.s100_pos-p.s0_pos)*x+p.s0_pos)-p.U_neg((p.s100_neg-p.s0_neg)*x+p.s0_neg); 
    soc_init = interp1(EMF,x,init_cond,'linear','extrap'); 
    if soc_init <0 || soc_init>1
        warning('Initial SOC not in the range of [0,1]')
    end
else
    soc_init = init_cond; 
end

if p.set_simp(6)
    cs(1:p.nn) = (cs100_neg-cs0_neg)*soc_init+cs0_neg; 
    cs(p.nn+1:p.nnp) = (cs100_pos-cs0_pos)*soc_init+cs0_pos; 
    cs(p.nnp+1:p.nnp+p.nn) = (cs100_neg-cs0_neg)*soc_init+cs0_neg; 
    cs(p.nnp+p.nn+1:2*p.nnp) = (cs100_pos-cs0_pos)*soc_init+cs0_pos; 
    cs = cs'; 
    cs_bar= p.Acs_bar*cs;
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
soc = soc_init; 
end

function [ p ] = fcn_system_matrices( p )
%-------------------------------------------------------------------------%
%- This function defines the system matrices that do not change in the ---%
%- Inner loop. The matrices that do change in the inner loop are placed --%
%- in their respective functions. ----------------------------------------%
%-------------------------------------------------------------------------%
if p.set_simp(6)
    p.Acs_hat = [-eye(p.nnp) zeros(p.nnp); -eye(p.nnp) eye(p.nnp)]; 
    p.Acs_hat_inv = inv(p.Acs_hat); 
else
Acs_jn = fcn_Acs_j(p.nrn);
p.Acs_jp = fcn_Acs_j(p.nrp);
Acs_n = sparse(kron(eye(p.nn),Acs_jn)*diag((1/p.dr_n^2)*p.Ds_neg*ones(p.nn,1)'*kron(eye(p.nn), ones(1,p.nrn))));
p.Acs_hat_n = p.dt*Acs_n-speye(p.nrn*p.nn);
p.Acs_hat_inv_n = inv(p.Acs_hat_n); 

Bcs_jn = -2*(p.dr_n+p.R_neg)/(p.dr_n*p.R_neg)*[zeros(p.nrn-1,1); 1];
Bcs_jp = -2*(p.dr_p+p.R_pos)/(p.dr_p*p.R_pos)*[zeros(p.nrp-1,1); 1];
Bcs_n = sparse(kron(eye(p.nn), Bcs_jn));
Bcs_p = sparse(kron(eye(p.np), Bcs_jp));
p.Bcs = blkdiag(Bcs_n,Bcs_p);
p.Bcs_hat = p.Bcs*p.dt; 
end

p.Bce = (diag(p.a_s.*(1-p.t_plus)./(p.eps_e(p.elec_range)))); 
p.Bce = sparse([p.Bce(1:p.nn,:); zeros(p.ns,p.nnp); p.Bce(p.nn+1:end,:)]);
p.Bce_hat = p.dt*p.Bce; 

Aphis_neg = fcn_Aw(p.sigma_eff(1:p.nn),ones(p.nn,1),p.dx(1:p.nn),p);
Aphis_pos = fcn_Aw(p.sigma_eff(p.nn+1:end),ones(p.np,1),p.dx(p.nns+1:end),p);
p.Aphis = sparse(blkdiag(Aphis_neg, Aphis_pos)); 

p.Bphis = sparse(diag(-p.a_s*p.F.*p.dx(p.elec_range))); 
p.Cphis = sparse([-1/p.A_surf; zeros(p.nnp-2,1); 1/p.A_surf]); 

p.Bphie = (diag(p.a_s*p.F.*p.dx(p.elec_range))); 
p.Bphie = sparse([p.Bphie(1:p.nn,:); zeros(p.ns,p.nnp); p.Bphie(p.nn+1:end,:)]);
p.Bphie(end) = 0; 

if p.set_simp(6)
    p.Acs_bar = [zeros(p.nnp) eye(p.nnp)]; 
else
Acs_temp_n = kron(eye(p.nn),[zeros(1,p.nrn-1) 1]); 
Acs_temp_p = kron(eye(p.np),[zeros(1,p.nrp-1) 1]);
p.Acs_bar = sparse(blkdiag(Acs_temp_n,Acs_temp_p)); 
end
Ace_temp = blkdiag(eye(p.nn), zeros(p.ns), eye(p.np)); 
p.Ace_bar = sparse([Ace_temp(1:p.nn,:); Ace_temp(p.nns+1:end,:)]); 
p.Aphie_bar = p.Ace_bar;

p.Bphis_inv = inv(p.Bphis); 

if p.set_simp(2)==2 || p.set_simp(2)==0
    p = De_matrices(p.ce0,p.T_amb,p); 
end

if p.set_simp(4)==2 || p.set_simp(4)== 0
    p = Ds_matrices(ones(p.np,1),p.T_amb,p); 
end

if p.set_simp(1)==2 || p.set_simp(1) == 0 
    p = kappa_matrices(p.ce0,p.T_amb,p); 
    p = nu_matrices(p.ce0,p.T_amb,p); 
end
end

function [ p ] = fcn_system_matrices2( p, cs_prevt, ce_prevt )
%-------------------------------------------------------------------------%
if p.set_simp(2)==1
    p = De_matrices(ce_prevt,p.T_amb,p); 
end

if p.set_simp(4)==1
    cs_bar = p.Acs_bar*cs_prevt; 
    p = Ds_matrices(cs_bar(p.nn+1:end)/p.cs_max_pos,p.T_amb,p); 
end

if p.set_simp(1)==1
    p = kappa_matrices(ce_prevt,p.T_amb,p); 
    p = nu_matrices(ce_prevt,p.T_amb,p); 
elseif p.set_simp(3)==1 && not(p.set_simp(1)==0)
    p = nu_matrices(ce_prevt,p.T_amb,p); 
end

if p.set_simp(4)==1 || p.set_simp(4)==2 || p.set_simp(4)==0
    if p.set_simp(6)==1
        p.Theta_cs = -p.Acs_hat_inv*[cs_prevt(1:p.nnp);zeros(p.nnp,1)]; 
    else
    p.Theta_cs = -p.Acs_hat_inv*cs_prevt; 
    end
    p.Theta_cs_bar = p.Acs_bar*p.Theta_cs; 
end
if p.set_simp(2)==1 || p.set_simp(2)==2 || p.set_simp(2)==0
    p.Theta_ce = -p.Ace_hat_inv*ce_prevt; 
    p.Theta_ce_bar = p.Ace_bar*p.Theta_ce;
end
end 

function p = De_matrices(ce,T,p)
    p.Ace = sparse(fcn_Aw(p.De_eff(ce,T),p.eps_e.*p.dx,p.dx,p));
    p.Ace_hat = (p.dt*p.Ace-speye(p.nx)); 
    p.Ace_hat_inv = inv(p.Ace_hat); 
    prem2 = p.Ace_hat_inv*(p.Bce_hat*p.Bphis_inv); 
    p.Gamma_ce = prem2*p.Cphis;
    p.Phi_ce = prem2*p.Aphis;
    p.Gamma_ce_bar = p.Ace_bar*p.Gamma_ce;
    p.Phi_ce_bar = p.Ace_bar*p.Phi_ce; 
end

function p = Ds_matrices(stoich,T,p)
    if p.set_simp(6)
        p.Ds = [p.Ds_neg*ones(p.nn,1); ones(p.np,1).*p.Ds_pos(stoich)]; 
        p.Bcs_hat = [-diag(3*p.dt./p.Rs); diag(p.Rs./(5*p.Ds))]; 
    else
        Acs_p = sparse(kron(eye(p.np),p.Acs_jp)*diag((1/p.dr_p^2)*(ones(p.np,1).*p.Ds_pos(stoich))'*kron(eye(p.np), ones(1,p.nrp))));
        Acs_hat_p = p.dt*Acs_p-speye(p.nrp*p.np);
        Acs_hat_inv_p = inv(Acs_hat_p); 
        p.Acs_hat_inv = blkdiag(p.Acs_hat_inv_n,Acs_hat_inv_p);  
    end
    prem1 = p.Acs_hat_inv*(p.Bcs_hat*p.Bphis_inv); 
    p.Gamma_cs = prem1*p.Cphis;
    p.Phi_cs = prem1*p.Aphis;
    p.Gamma_cs_bar = p.Acs_bar*p.Gamma_cs;
    p.Phi_cs_bar = p.Acs_bar*p.Phi_cs; 
end

function p = kappa_matrices(ce,T,p)
    Aphie = (fcn_Aw(p.kappa_eff(ce,T),ones(p.nx,1),p.dx,p));
    p.Aphie = Aphie;
    p.Aphie(end,end-1:end) = [0 1]; 
    p.Aphie_inv = inv(p.Aphie); 
    prem4 = p.Aphie_inv*(p.Bphie*p.Bphis_inv); 
    p.Gamma_phie = prem4*p.Cphis;
    p.Phi_phie= prem4*p.Aphis;
    p.Gamma_phie_bar = p.Aphie_bar*p.Gamma_phie;
    p.Phi_phie_bar = p.Aphie_bar*p.Phi_phie; 
end

function p = nu_matrices(ce,T,p)
    p.Dphie = fcn_Aw(p.nu(ce,T),ones(p.nx,1),p.dx,p);
    p.Dphie(end,end-1:end) = [0 0]; 
    p.Pi_phie= -p.Aphie_inv*p.Dphie; 
%     p.Pi_phie = zeros(p.nx,p.nx); 
    p.Pi_phie_bar = p.Aphie_bar*p.Pi_phie; 
    
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
function [ p ] = fcn_system_vectors( p,t)
%-------------------------------------------------------------------------%
%- This function converts the specified parameters into vectors. ---------%
%-------------------------------------------------------------------------%
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

p.k0 = [p.k0_neg*ones(p.nn,1); p.k0_pos*ones(p.np,1)];
p.Rs = [p.R_neg*ones(p.nn,1); p.R_pos*ones(p.np,1)]; 
p.R_f = [p.Rf_neg*ones(p.nn,1); p.Rf_pos*ones(p.np,1)];

p.eps_e = [p.epse_neg*ones(p.nn,1); p.epse_sep*ones(p.ns,1); p.epse_pos*ones(p.np,1)];
p.eps_s = [p.epss_neg*ones(p.nn,1); p.epss_pos*ones(p.np,1)];
p.dx = [p.dx_n*ones(p.nn,1); p.dx_s*ones(p.ns,1); p.dx_p*ones(p.np,1)]; 
p.brug = [p.p_neg*ones(p.nn,1); p.p_sep*ones(p.ns,1); p.p_pos*ones(p.np,1)]; 
p.stoich_min = [p.s100_neg*ones(p.nn*p.nrn,1); p.s100_pos*ones(p.np*p.nrp,1)]; 
p.stoich_max = [p.s100_neg*ones(p.nn*p.nrn,1); p.s0_pos*ones(p.np*p.nrp,1)]; 

p.sigma_eff_neg = p.sigma_neg.*p.epss_neg;                                    %Effective solid-phase conductivity at the neg. electrode [S/m]
p.sigma_eff_pos = p.sigma_pos*p.epss_pos;                                    %Effective solid-phase conductivity at the pos. electrode [S/m]   
p.sigma_eff = [p.sigma_eff_neg*ones(p.nn,1); p.sigma_eff_pos*ones(p.np,1)];
p.a_s_neg = 3*p.epss_neg/p.R_neg;                                              %Specific interfacial surface area at the neg. electrode [1/m]
p.a_s_pos = 3*p.epss_pos/p.R_pos;                                              %Specific interfacial surface area at the pos. electrode [1/m]
p.a_s = [p.a_s_neg*ones(p.nn,1); p.a_s_pos*ones(p.np,1)];

parnames = {'kappa' 'De' 'dlnfdx' 'Ds_pos'}; 
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

p.kappa_eff = @(c,T) p.kappa(c,T).*p.eps_e.^p.brug;     %[S/m] 
p.nu = @(c,T) 2*p.R*T*p.kappa_eff(c,T)*(p.t_plus-1)/p.F.*(1+p.dlnfdx(c,T)); 
p.De_eff = @(c,T) p.De(c,T).*p.eps_e.^p.brug; 

% p.elec_range specifies the range of the electrodes throughout the cell.
p.elec_range = [linspace(1,p.nn,p.nn) linspace(p.nns+1,p.nx, p.np)];

p.cs_bar_max = [p.cs_max_neg*ones(p.nn,1); p.cs_max_pos*ones(p.np,1)];
p.cs_max = [p.cs_max_neg*ones(p.nn*p.nrn,1); p.cs_max_pos*ones(p.np*p.nrp,1)];
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
function print_loop(t1,time_max,i_app,k,sim_time)
dispstat(sprintf(['-------------------------------------------\nTime: %g/',num2str(time_max),...
    ' (%.0f %%) \ni_app: %g A\nNumber of iterations for convergence: %d \nTime elapsed: ',num2str(sim_time,'%2.2f'),'\n-------------------------------------------'],...
    t1,t1/time_max*100,i_app,k));
end