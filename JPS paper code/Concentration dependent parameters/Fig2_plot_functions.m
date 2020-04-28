clear all; close all
fontsize = 13; 
load('Ds_Safari_Discharge.mat')
load Ds_Farkhondeh.mat

figure(1)
subplot(2,2,4)
plot(stoich_C,log10(Ds_C),'k','LineWidth',2)
hold on
plot(stoich,log10(Ds),'r','LineWidth',2)
Ds_fit = 10.^(-20.26+534.9*(stoich-0.5).^8+2.263*(stoich-0.5).^2); 
plot(stoich,log10(Ds_fit),'r--','LineWidth',2)

grid on
xlim([0 1]); 
ylim([-20.4 -17])
set(gca,'ytick',-20:1:17)
ax = gca; 
ax.YGrid = 'on'; 
ax.XGrid = 'on';
xlabel('$c_s/c_{s,\rm{max}} \ \rm{[-]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
ylabel('$\log_{10}(D_{s,p}) \ \rm{[m^2/s]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
lgd1 = legend({'$D_{s,p}$ from [35]','$D_{s,p}$ from [24]', 'Approx. of [24]'},'Interpreter','latex','FontWeight','bold');


subplot(2,2,3)
t_plus = 0.363; 
nu = @(c,T) (0.601-0.24*(c/1000).^0.5+0.983.*(1-0.0052*(T-294))*(c/1000).^1.5); 
T = 293; 
ce = 0:1:2000; 
plot(ce,nu(ce,T),'k','LineWidth',2)
xlabel('$c_e \ \rm{[mol/m^3]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
ylabel('$\nu \ \rm[-]$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
grid on
lgd2 = legend({'$\nu$ from [21]'},'Interpreter','latex','FontWeight','bold');

subplot(2,2,1)
kappa = @(c,T) 0.1*(c/1000).*(-10.5+0.668*(c/1000)+0.494*(c/1000).^2+0.074*T-0.0178*(c/1000)*T-8.86*10^(-4)*(c/1000).^2*T-6.96*10^(-5)*T^2+2.81*10^(-5)*(c/1000)*T^2).^2; 
kappa2 = @(c) 15.8e-4*c.*exp(0.85*(1e-3*c).^1.4); 
kappa3 = @(c_e) 0.0911+1.9101*c_e/1e3 - 1.052*(c_e/1e3).^2 + 0.1554*(c_e/1e3).^3;
kappa4 = @(c) 15.8e-4*c.*exp(-0.85*(1e-3*c).^1.4); 
plot(ce,kappa(ce,T),'k','LineWidth',2)
hold on
plot(ce,kappa3(ce),'b','LineWidth',2)
hold on
plot(ce,kappa4(ce),'r','LineWidth',2)
plot(ce,kappa2(ce)/10,'r--','LineWidth',2)
xlabel('$c_e \ \rm{[mol/m^3]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
ylabel('$\kappa \ \rm{[S/m]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
grid on
lgd3 =legend({'$\kappa$ from [21]','$\kappa$ from [34]','$\kappa$ from [38]','$\kappa/10$ from [31]'},'Interpreter','latex','FontWeight','bold');

subplot(2,2,2)
De = @(ce,T) 10e-4*10.^(-4.43-54./(T-229-5*(ce/1000))-0.22*(ce/1000)); 
De2 = @(c_e) 5.34e-10*exp(-0.65*c_e/1e3); 
plot(ce,De(ce,T),'k','LineWidth',2)
hold on
plot(ce,10*De2(ce),'r','LineWidth',2)
xlabel('$c_e \ \rm{[mol/m^3]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
ylabel('$D_e \ \rm{[m^2/s]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
grid on
lgd4=legend({'$D_e$ from [21]','$10D_e$ from [34]'},'Interpreter','latex','FontWeight','bold');

set(gcf, 'Position',  [20, 20, 850, 750])
set(findall(gcf,'-property','FontSize'),'FontSize',20)

lgd1.FontSize = 15; 
lgd2.FontSize = 16; 
lgd3.FontSize = 16;
lgd4.FontSize = 16; 