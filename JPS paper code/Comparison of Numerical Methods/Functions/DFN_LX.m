%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to simulate the DFN model in which the numerical method 
% presented in [1] is implemented
%
% Model Simplifications and Its Impact on Computational Complexity for an 
% Electrochemistry-Based Battery Modeling Toolbox
%
% Authors: Z. Khalik, M.C.F. Donkers, H.J. Bergveld
%
% This file is licensed under the BSD 3-Clause License
%
% References
% [1] Xia et al, A computationally efficient implementation of a full and
% reduced-order electrochemistry-based model for Li-ion batters, Applied
% Energy, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function out = DFN_LX(input_current,tf,soc_init,param)
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
p = fcn_system_vectors(p);
p = fcn_system_matrices(p);
p = fcn_bv_derivatives(p);

[cs(:,1), ce(:,1), phis(:,1), phie(:,1)] = init(p,soc_init);

%- prev indicates the condition of the state at k-1
cs_prev = cs(:,1);
ce_prev = ce(:,1);
phis_prev = phis(:,1);
phie_prev = phie(:,1);
cs_prevt = cs_prev;
ce_prevt = ce_prev;
soc_prevt = soc_init;

p.gamma_init = 1;
gamma = p.gamma_init;
                                                  %Measure simulation time for benchmark purposes
dispstat('Starting simulation...','keepthis','timestamp')
warning_set = 0;
t_vec = p.dt; 
end_simulation = 0; 
solution_time = 0;
not_converged = 0;

%% Simulation
for t=1:1e6
    t1=t*p.dt;
    tic()
    %- Modified Newton's Method ------------------------------------------%
   for k=1:p.iter_max
       [F_cs, J_cs] = fcn_cs(cs_prev,ce_prev,phis_prev,phie_prev,cs_prevt,p);
       cs(:,t) = cs_prev-gamma*(J_cs\F_cs); 
       conv_check(1) = norm((cs(:,t)-cs_prev),2);
       
       [F_ce, J_ce] = fcn_ce(cs(:,t),ce_prev,phis_prev,phie_prev,ce_prevt,p);
       ce(:,t) = ce_prev-gamma*(J_ce\F_ce); 
       conv_check(2) = norm((ce(:,t)-ce_prev),2);
       
       [F_phie, J_phie] = fcn_phie(cs(:,t),ce(:,t),phis_prev,phie_prev,p);
       phie(:,t) = phie_prev-gamma*(J_phie\F_phie);
       conv_check(3) = norm((phie(:,t)-phie_prev),2);
       
       [F_phis, J_phis] = fcn_phis(cs(:,t),ce(:,t),phis_prev,phie(:,t),i_app(t),p);
       phis(:,t) = phis_prev-gamma*(J_phis\F_phis); 
       conv_check(4) = norm((phis(:,t)-phis_prev),2);

       cs_prev = cs(:,t);
       ce_prev = ce(:,t);
       phis_prev = phis(:,t);
       phie_prev = phie(:,t);
       
       % If algorithm doesn't converge, then display a warning
       if k==p.iter_max
           warning('NOT CONVERGED at time %.0f', t1) 
           end_simulation=1;
       end    
       
       if(all(conv_check(2)' <=p.tol))
           gamma = p.gamma;
           print_loop(t1,time_max,i_app(t),k,solution_time)
           break
       end
   end  
   % Abort simulation if algorithm is not converged
   [jn(:,t), ~,U(:,t),eta(:,t)] = fcn_bv(cs_prev,ce_prev,phis_prev,phie_prev, p ); 
   solution_time = solution_time+toc();
    if not_converged ==1
        break
    end
    V(t) = phis(end,t)+i_app(t)/p.A_surf*p.dx(end)*0.5/p.sigma_eff(end)...
    -phis(1,t)-i_app(t)/p.A_surf*p.dx(1)*0.5/p.sigma_eff(1)...
    +(p.R_cc/p.A_surf)*i_app(t);
    T(t) = p.T_amb; 
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
if warning_set==0
    if p.verbose==0
    dispstat(sprintf('Finished the simulation in %2.2f s \n',out.sim_time),'keepthis','timestamp');
    end
else
    dispstat(sprintf('Finished the simulation in %2.2f s with warning %d \n',out.sim_time,warning_set),'keepthis','timestamp');
end
i_app = i_app(1:end-1);
n_t = t; 
p.L = p.delta_neg+p.delta_sep+p.delta_pos; 
% Store states
if not(p.output_interpolated_states)
out.x = [p.dx_n/2*(1:2:(p.nn*2-1)) p.delta_neg+[p.dx_s/2*(1:2:(p.ns*2-1))]...
        p.delta_neg+p.delta_sep+p.dx_p/2*(1:2:(p.np*2-1))]/p.L; 
out.phis = [phis(1:p.nn,:); NaN(p.ns,n_t); phis(p.nn+1:end,:)]'; 
out.ce = ce'; 
out.cs = cs';  
out.phie = phie'; 
if p.set_simp(6)
    out.cs_bar = [cs(p.nnp+(1:p.nn),:); NaN(p.ns,n_t); cs(p.nnp+(p.nn+1:p.nnp),:)]';
else
    out.cs_bar = [cs(p.nrn:p.nrn:p.nrn*p.nn,:); NaN(p.ns,n_t); cs(p.nrn*p.nn+p.nrp:p.nrp:end,:)]';
end
out.stoich = [out.cs_bar(:,1:p.nn)/p.cs_max_neg NaN(n_t,p.ns) out.cs_bar(:,p.nns+1:end)/p.cs_max_pos]; 
out.jn = [jn(1:p.nn,:); NaN(p.ns,n_t); jn(p.nn+1:end,:)]';  
out.U = [U(1:p.nn,:); NaN(p.ns,n_t); U(p.nn+1:end,:)]';  
out.eta = [eta(1:p.nn,:); NaN(p.ns,n_t); eta(p.nn+1:end,:)]'; 
out.V = V'; 
out.i_app = i_app';
out.T = T'; 
out.p = p; 
out.t = t_vec(1:end-1)'; 
out.soc = soc';

else
% Store extrapolated states if required
x_n = [0 p.dx_n/2*(1:2:(p.nn*2-1)) p.delta_neg]; 
x_s = p.delta_neg+[p.dx_s/2*(1:2:(p.ns*2-1))]; 
x_p = p.delta_neg+p.delta_sep+[0 p.dx_p/2*(1:2:(p.np*2-1)) p.delta_pos]; 
x_s_ip = [x_n x_s x_p]'; 
x_e_ip = [x_n(1:end-1) x_s x_p(2:end)]; 
phis_neg = [phis(1,:)-i_app*p.dx(1)*0.5/p.sigma_eff(1); phis(1:p.nn,:); phis(p.nn,:)]; 
phis_sep = NaN(p.ns,size(phis_neg,2)); 
phis_pos = [phis(p.nn+1,:); phis(p.nn+1:end,:); phis(p.nnp,:)+i_app*p.dx(end)*0.5/p.sigma_eff(p.nnp)]; 
phis_ip = [phis_neg; phis_sep; phis_pos]; 
phie_ip = [phie(1,:); phie(1:p.nn,:);interp1(x_e_ip(p.nn+1:p.nn+2),phie(p.nn:p.nn+1,:),p.delta_neg); phie(p.nn+1:p.nns,:);...
        interp1(x_e_ip(p.nns+1:p.nns+2),phie(p.nns:p.nns+1,:),p.delta_neg+p.delta_sep);phie(p.nns+1:end,:);phie(end,:)];
ce_ip = [ce(1,:); ce(1:p.nn,:);interp1(x_e_ip(p.nn+1:p.nn+2),ce(p.nn:p.nn+1,:),p.delta_neg); ce(p.nn+1:p.nns,:);...
        interp1(x_e_ip(p.nns+1:p.nns+2),ce(p.nns:p.nns+1,:),p.delta_neg+p.delta_sep);ce(p.nns+1:end,:);ce(end,:)];
if p.set_simp(6)
    cs_bar_neg = [NaN(1,t); cs(p.nnp+1:p.nnp+p.nn,:); NaN(1,t)]; 
    cs_bar_sep = NaN(p.ns,t); 
    cs_bar_pos = [NaN(1,t); cs(p.nnp+p.nn+1:end,2:t); NaN(1,t)]; 
else
    cs_bar_neg = [NaN(1,t); cs(p.nrn:p.nrn:p.nrn*p.nn,:); NaN(1,t)]; 
    cs_bar_sep = NaN(p.ns,t); 
    cs_bar_pos = [NaN(1,t); cs(p.nrn*p.nn+p.nrp:p.nrp:end,:); NaN(1,t)]; 
end
cs_bar_ip = [cs_bar_neg; cs_bar_sep; cs_bar_pos];
stoich_ip = [cs_bar_neg/p.cs_max_neg; cs_bar_sep; cs_bar_pos/p.cs_max_pos]; 
U_neg = [NaN(1,t); U(1:p.nn,1:t); NaN(1,t)]; 
U_sep = NaN(p.ns,t); 
U_pos = [NaN(1,t); U(p.nn+1:end,1:t); NaN(1,t)]; 
U_ip = [U_neg; U_sep; U_pos]; 
jn_neg = [NaN(1,t); jn(1:p.nn,1:t); NaN(1,t)]; 
jn_sep = NaN(p.ns,t); 
jn_pos = [NaN(1,t); jn(p.nn+1:end,1:t); NaN(1,t)]; 
jn_ip = [jn_neg; jn_sep; jn_pos]; 
eta_neg = [NaN(1,t); eta(1:p.nn,1:t); NaN(1,t)]; 
eta_sep = NaN(p.ns,t); 
eta_pos = [NaN(1,t); eta(p.nn+1:end,1:t); NaN(1,t)]; 
eta_ip = [eta_neg; eta_sep; eta_pos]; 

out.x = x_s_ip'/p.L; 
out.phis = phis_ip'; 
out.ce = ce_ip';
out.cs = cs'; 
out.phie = phie_ip'; 
out.cs_bar = cs_bar_ip';
out.stoich = stoich_ip'; 
out.jn = jn_ip'; 
out.U = U_ip'; 
out.eta = eta_ip'; 
out.V = V'; 
out.i_app = i_app';
out.T = T'; 
out.p = p; 
out.t = t_vec(1:end-1)'; 
out.soc = soc';
end
end

%% Functions
%-------------------------------------------------------------------------%
%-- Functions for the DFN model ------------------------------------------%
%-------------------------------------------------------------------------%
function [ jn, i0, U,eta] = fcn_bv(cs,ce,phis,phie, p )
%-------------------------------------------------------------------------%
%- Compute j_li using a Butler-Volmer equation ---------------------------%
%- See Eq. (6) in Xia et al. (2017) --------------------------------------%
%-------------------------------------------------------------------------%
%- Concentration of the solid phase at r = R_s

cs_bar = p.Acs_bar*cs; 

%- Exchange current density, Eq. (7) in Xia et al. (2017)
i0 = p.k0.*ce(p.elec_range).^p.alpha_a.*(p.cs_bar_max-cs_bar).^p.alpha_a.*cs_bar.^p.alpha_c; 

%- Open circuit potential (see Table 1 in Xia et al. (2017))
w = cs_bar(1:p.nn)/p.cs_max_neg; 
z = cs_bar(p.nn+1:end)/p.cs_max_pos; 
U = [p.U_neg(w); p.U_pos(z)]; 
%- Over potential, Eq. (8)

% Linearized jn, also acts as first guess for nonlinearized jn
jn = (phis-phie(p.elec_range)-U)./(p.R*p.T./i0+p.F*p.R_f); 
eta = phis-phie(p.elec_range)-U-p.R_f*p.F.*jn; 
    
if p.set_simp(5)==0
    exp1 = exp((p.alpha_a*p.F*eta)/(p.R*p.T));
    exp2 = exp(-(p.alpha_c*p.F*eta)/(p.R*p.T));
    if any(p.R_f>0)
% Solve jn using Newton's method
        converged = 0;
        while(not(converged))
            jn_prev = jn; 
            F_jn = jn_prev-(i0/p.F).*(exp1-exp2); 
            J_jn = eye(p.nnp)+diag((i0.*p.R_f*p.F/(p.R*p.T)).*(p.alpha_a*exp1+p.alpha_c*exp2)); 
            jn = jn_prev-(J_jn\F_jn); 
            if norm((jn-jn_prev),2)<1e-10
                converged = 1;
            end
        end
    else
        jn = i0/p.F.*(exp1-exp2); 
    end
end

end
function [ F_ce, J_ce ] = fcn_ce(cs,ce,phis,phie,ce_prevt,p)
%-------------------------------------------------------------------------%
%- Compute F_ce and J_ce, Eq. (14c) --------------------------------------%
%-------------------------------------------------------------------------%
[jn,i0,U] = fcn_bv(cs,ce,phis,phie,p);
dj_dce = p.dj_dce(ce,phis,phie,i0,U);

if p.set_simp(2)==0
    p.Ace = sparse(fcn_Aw(p.De_eff(ce,p.T),p.eps_e.*p.dx,p.dx,p));
    p.Ace_hat = (p.dt*p.Ace-speye(p.nx)); 
end

F_ce = p.Ace_hat*ce+p.Bce_hat*jn+ce_prevt;
J_ce = p.Ace_hat+spdiags(p.Bce_hat*dj_dce,0,p.nx,p.nx);
end
function [ F_cs, J_cs ] = fcn_cs(cs,ce,phis,phie,cs_prevt,p)
%-------------------------------------------------------------------------%
%- Compute F_cs and J_cs, Eq. (14b) --------------------------------------%
%-------------------------------------------------------------------------%
[jn ,i0,U] = fcn_bv(cs,ce,phis,phie,p);
cs_bar = p.Acs_bar*cs; 
dj_dcs_bar = p.dj_dcs_bar(cs_bar,phis,phie,jn,i0,U); 

if p.set_simp(6)==1
    F_cs = p.Acs_hat*cs+p.Bcs_hat*jn-cs_prevt;
    J_cs = p.Acs_hat+spdiags(p.Bcs_hat*dj_dcs_bar,0,2*p.nn+2*p.np,2*p.nn+2*p.np);
else
    F_cs = p.Acs_hat*cs+p.Bcs_hat*jn+cs_prevt;
    J_cs = p.Acs_hat+spdiags(p.Bcs_hat*dj_dcs_bar,0,p.nrn*p.nn+p.nrp*p.np,p.nrn*p.nn+p.nrp*p.np);
end
end
function [ F_phie, J_phie ] = fcn_phie(cs,ce,phis,phie,p)
%-------------------------------------------------------------------------%
%- Compute F_phie and J_phie, Eq. (14d) ----------------------------------%
%-------------------------------------------------------------------------%
[jn, i0] = fcn_bv(cs,ce,phis,phie,p);
dj_dphie = p.dj_dphie(i0);

if p.set_simp(1)==0
    p.Dphie = fcn_Aw(p.nu(ce,p.T),ones(p.nx,1),p.dx,p);
    p.Dphie(end,end-1:end) = [0 0]; 
    p.Aphie = fcn_Aw(p.kappa_eff(ce,p.T),ones(p.nx,1),p.dx,p);
    p.Aphie(end,end-1:end) = [0 1]; 
elseif p.set_simp(3)==0
    p.Dphie = fcn_Aw(p.nu(ce,p.T),ones(p.nx,1),p.dx,p);
    p.Dphie(end,end-1:end) = [0 0]; 
end


F_phie = p.Aphie*phie+p.Dphie*log(ce)+p.Bphie*jn;
J_phie = p.Aphie+spdiags(p.Bphie*dj_dphie,0,p.nx,p.nx);
end
function [ F_phis, J_phis ] = fcn_phis(cs,ce,phis,phie, i_app,p)
%-------------------------------------------------------------------------%
%- Compute F_phis and J_phis, Eq. (14e) ----------------------------------%
%-------------------------------------------------------------------------%
[jn, i0] = fcn_bv(cs,ce,phis,phie,p);
dj_dphis = p.dj_dphis(i0);

F_phis = p.Aphis*phis+p.Bphis*jn+p.Cphis*i_app;

J_phis = p.Aphis+spdiags(p.Bphis*dj_dphis,0,p.nnp,p.nnp);
end
%Other functions

function p = fcn_bv_derivatives(p)
p.di0_dcs_bar = @(cs_bar, i0) (-p.alpha_a./(p.cs_bar_max-cs_bar)+p.alpha_c./cs_bar).*i0;
p.dU = @(cs_bar) [p.dU_neg(cs_bar(1:p.nn)/p.cs_max_neg)/p.cs_max_neg; p.dU_pos(cs_bar(p.nn+1:end)/p.cs_max_pos)/p.cs_max_pos]; 

p.dj_dcs_bar=@(cs_bar,phis,phie,jn,i0,U) -p.dU(cs_bar)./(p.R*p.T./i0+p.F*p.R_f)+...
    (phis-phie(p.elec_range)-U)*p.R*p.T.*p.di0_dcs_bar(cs_bar,i0)./(p.R*p.T+p.F*p.R_f.*i0).^2; 

p.di0_dce = @(ce,i0) p.alpha_a./ce(p.elec_range); 

p.dj_dce = @(ce,phis,phie,i0,U) (phis-phie(p.elec_range)-U)*p.R*p.T.*p.di0_dce(ce,i0)./(p.R*p.T+p.F*p.R_f.*i0).^2; 

p.dj_dphie = @(i0) -1./(p.R*p.T./i0+p.F*p.R_f); 

p.dj_dphis = @(i0) 1./(p.R*p.T./i0+p.F*p.R_f); 

end

function [ cs, ce, phis, phie ] = init(p,soc_init)
%-------------------------------------------------------------------------%
%- Initial conditions for the states cs, ce, phis, phie. -----------------%
%-------------------------------------------------------------------------%
cs0_neg = p.s0_neg*p.cs_max_neg;
cs100_neg = p.s100_neg*p.cs_max_neg;
cs0_pos = p.s0_pos*p.cs_max_pos;
cs100_pos = p.s100_pos*p.cs_max_pos;

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

phie = zeros(p.nx,1); 
end

function [ A1_j ] = fcn_Acs_j(nr )
temp1 = linspace(0,nr-2,nr-1)./linspace(1,nr-1,nr-1);
temp2 = [2 linspace(2,nr-1,nr-2)./linspace(1,nr-2,nr-2)];

temp3 = diag(temp1,-1)-2*eye(nr)+diag(temp2,1);
temp3(nr,nr-1) = 2;
A1_j = temp3;
end
function [ Aw] = fcn_Aw(alpha1, alpha0, dx,p)
switch p.fvm_method
    case 1 %Harmonic mean
        alpha1_face= alpha1(2:end).*alpha1(1:end-1).*(dx(2:end)+dx(1:end-1))./(alpha1(2:end).*dx(1:end-1)+alpha1(1:end-1).*dx(2:end)); 
    case 2 %Linear variation
        alpha1_face = (alpha1(2:end).*dx(1:end-1)+alpha1(1:end-1).*dx(2:end))./(dx(2:end)+dx(1:end-1)); 
    case 3 %Mean
        alpha1_face = (alpha1(2:end).*dx(2:end)+alpha1(1:end-1).*dx(1:end-1))./(dx(2:end)+dx(1:end-1)); 
end
gamma = alpha1_face./(dx(2:end)+dx(1:end-1)); 

Am = diag(gamma,-1) + diag([-gamma(1); -gamma(2:end)-gamma(1:end-1); -gamma(end)])+diag(gamma,1); 
Aw0 = diag(2./(alpha0)); 
Aw = Aw0*Am; 
end
function [ p ] = fcn_system_matrices( p )
%-------------------------------------------------------------------------%
%- This function defines the system matrices that do not change in the ---%
%- Inner loop. The matrices that do change in the inner loop are placed --%
%- in their respective functions. ----------------------------------------%
%-------------------------------------------------------------------------%
if p.set_simp(6)
    p.Acs_bar = [zeros(p.nnp) eye(p.nnp)]; 
else
Acs_temp_n = kron(eye(p.nn),[zeros(1,p.nrn-1) 1]); 
Acs_temp_p = kron(eye(p.np),[zeros(1,p.nrp-1) 1]);
p.Acs_bar = blkdiag(Acs_temp_n,Acs_temp_p); 
end

if p.set_simp(6)
    p.Acs_hat = [eye(p.nnp) zeros(p.nnp); -eye(p.nnp) eye(p.nnp)]; 
    p.Acs_hat_inv = inv(p.Acs_hat); 
    p.Ds = [p.Ds_neg*ones(p.nn,1); ones(p.np,1).*p.Ds_pos((p.s100_pos+p.s0_pos)/2)]; 
    p.Bcs_hat = [diag(3*p.dt./p.Rs); diag(p.Rs./(5*p.Ds))]; 
else
Acs_jn = fcn_Acs_j(p.nrn);
p.Acs_jp = fcn_Acs_j(p.nrp);
Acs_n = sparse(kron(eye(p.nn),Acs_jn)*diag((1/p.dr_n^2)*p.Ds_neg*ones(p.nn,1)'*kron(eye(p.nn), ones(1,p.nrn))));
p.Acs_hat_n = p.dt*Acs_n-speye(p.nrn*p.nn);
Acs_p = sparse(kron(eye(p.np),p.Acs_jp)*diag((1/p.dr_p^2)*(ones(p.np,1).*p.Ds_pos((p.s100_pos+p.s0_pos)/2))'*kron(eye(p.np), ones(1,p.nrp))));
Acs_hat_p = p.dt*Acs_p-speye(p.nrp*p.np);
p.Acs_hat = blkdiag(p.Acs_hat_n,Acs_hat_p); 

Bcs_jn = -2*(p.dr_n+p.R_neg)/(p.dr_n*p.R_neg)*[zeros(p.nrn-1,1); 1];
Bcs_jp = -2*(p.dr_p+p.R_pos)/(p.dr_p*p.R_pos)*[zeros(p.nrp-1,1); 1];
Bcs_n = sparse(kron(eye(p.nn), Bcs_jn));
Bcs_p = sparse(kron(eye(p.np), Bcs_jp));
p.Bcs = blkdiag(Bcs_n,Bcs_p);
p.Bcs_hat = p.Bcs*p.dt; 
end

if p.set_simp(2)==2
    p.Ace = sparse(fcn_Aw(p.De_eff(p.ce0,p.T),p.eps_e.*p.dx,p.dx,p));
    p.Ace_hat = (p.dt*p.Ace-speye(p.nx)); 
end

p.Bce = (diag(p.a_s.*(1-p.t_plus)./(p.eps_e(p.elec_range)))); 
p.Bce = sparse([p.Bce(1:p.nn,:); zeros(p.ns,p.nnp); p.Bce(p.nn+1:end,:)]);
p.Bce_hat = p.dt*p.Bce; 

Aphis_neg = fcn_Aw(p.sigma_eff(1:p.nn),ones(p.nn,1),p.dx(1:p.nn),p);
Aphis_pos = fcn_Aw(p.sigma_eff(p.nn+1:end),ones(p.np,1),p.dx(p.nns+1:end),p);
p.Aphis = sparse(blkdiag(Aphis_neg, Aphis_pos)); 

p.Bphis = sparse(diag(-p.a_s*p.F.*p.dx(p.elec_range))); 
p.Cphis = sparse([-1/p.A_surf; zeros(p.nnp-2,1); 1/p.A_surf]); 

if p.set_simp(3)==2
    p.Dphie = fcn_Aw(p.nu(p.ce0,p.T),ones(p.nx,1),p.dx,p);
    p.Dphie(end,end-1:end) = [0 0]; 
    p.Aphie = fcn_Aw(p.kappa_eff(p.ce0,p.T),ones(p.nx,1),p.dx,p);
    p.Aphie(end,end-1:end) = [0 1]; 
elseif p.set_simp(1)==2
    p.Aphie = fcn_Aw(p.kappa_eff(p.ce0,p.T),ones(p.nx,1),p.dx,p);
    p.Aphie(end,end-1:end) = [0 1]; 
end

p.Bphie = (diag(p.a_s*p.F.*p.dx(p.elec_range))); 
p.Bphie = sparse([p.Bphie(1:p.nn,:); zeros(p.ns,p.nnp); p.Bphie(p.nn+1:end,:)]);
p.Bphie(end) = 0; 
end
function [ p ] = fcn_system_vectors( p)
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
p.dx = [p.dx_n*ones(p.nn,1); p.dx_s*ones(p.ns,1); p.dx_p*ones(p.np,1)]; 
p.brug = [p.p_neg*ones(p.nn,1); p.p_sep*ones(p.ns,1); p.p_pos*ones(p.np,1)]; 

p.sigma_eff_neg = p.sigma_neg*p.epss_neg;                                    %Effective solid-phase conductivity at the neg. electrode [S/m]
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
                p.(parnames{i})= p.(parnames{i})((p.s100_pos+p.s0_pos)/2,p.T);
            else
                p.(parnames{i}) = p.(parnames{i})(p.ce0,p.T); 
            end
        end
        p.(parnames{i}) = @(c,T) p.(parnames{i}); 
        p.set_simp(i) = 2; 
    end
end

% if not(isa(p.Ds_pos,'function_handle')) || p.set_simp(4)==2
%     if isa(p.Ds_pos,'function_handle')
%         p.Ds_pos= p.Ds_pos((p.s100_pos+p.s0_pos)/2,p.T); 
%     end
%     p.Ds_pos = @(stoich,T) ones(length(stoich),1)*p.Ds_pos; 
%     p.set_simp(4) = 2; 
% end

p.kappa_eff = @(c,T) p.kappa(c,T).*p.eps_e.^p.brug;     %[S/m] 
p.nu = @(c,T) 2*p.R*T*p.kappa_eff(c,T)*(p.t_plus-1)/p.F.*(1+p.dlnfdx(c,T)); 
p.De_eff = @(c,T) p.De(c,T).*p.eps_e.^p.brug; 

% p.elec_range specifies the range of the electrodes throughout the cell.
p.elec_range = [linspace(1,p.nn,p.nn) linspace(p.nns+1,p.nx, p.np)];

p.cs_bar_max = [p.cs_max_neg*ones(p.nn,1); p.cs_max_pos*ones(p.np,1)];
p.cs_max = [p.cs_max_neg*ones(p.nn*p.nrn,1); p.cs_max_pos*ones(p.np*p.nrp,1)];
p.stoich_min = [p.s0_neg*ones(p.nn*p.nrn,1); p.s100_pos*ones(p.np*p.nrp,1)]; 
p.stoich_max = [p.s100_neg*ones(p.nn*p.nrn,1); p.s0_pos*ones(p.np*p.nrp,1)]; 
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
