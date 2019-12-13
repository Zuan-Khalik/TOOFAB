function make_plots(out,eval_time,plotcolor,plotcolor_V,fig_handle,set_scaleheight)
figure(fig_handle)
fontsize = 13; 
scaleheight = 0.6; 
subplot(3,2,1:2)
plot(out.t,out.V,plotcolor_V,'LineWidth',2)
hold on
grid on
xlabel('Time [s]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
ylabel('$V_t \ \mathrm{[V]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
xlim([0 out.t(end)])

subplot(3,2,3)
plot(out.x,out.phie(eval_time,:),plotcolor,'LineWidth',2); 
grid on
hold on
ylabel('$\phi_e \ \mathrm{[V]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
set(gca,'xtick',0:0.2:1)
if set_scaleheight
A = get(gca,'position');          % gca points at the second one
A(1,4) = A(1,4) *scaleheight;              % reduce the height by half
A(1,2) = A(1,2) +A(1,4)/scaleheight*(1-scaleheight);         % change the vertical position
set(gca,'position',A);     
end

subplot(3,2,4)
plot(out.x,out.ce(eval_time,:)/out.p.ce0,plotcolor,'LineWidth',2)
grid on
hold on
ylabel('$c_e/c_{e,0} \ \mathrm{[-]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
set(gca,'xtick',0:0.2:1)
if set_scaleheight
A = get(gca,'position');          % gca points at the second one
A(1,4) = A(1,4) *scaleheight;              % reduce the height by half
A(1,2) = A(1,2) +A(1,4)/scaleheight*(1-scaleheight);         % change the vertical position
set(gca,'position',A);     
end

subplot(3,2,5)
plot(out.x,out.stoich(eval_time,:),plotcolor,'LineWidth',2)
grid on
hold on
xlabel('x/L','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
ylabel('$\bar{c}_s/c_{s,\mathrm{max}} \ \mathrm{[-]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
set(gca,'xtick',0:0.2:1)
if set_scaleheight
A = get(gca,'position');          % gca points at the second one
A(1,4) = A(1,4) *scaleheight;              % reduce the height by half
A(1,2) = A(1,2) +2*A(1,4)/scaleheight*(1-scaleheight);         % change the vertical position
set(gca,'position',A);     
end

% eta = out.phis-out.phie-out.U;
eta = out.eta; 
subplot(3,2,6)
plot(out.x,eta(eval_time,:),plotcolor,'LineWidth',2)
grid on
hold on
xlim([0 1])
xlabel('x/L','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
ylabel('$\eta \ [V]$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
set(gca,'xtick',0:0.2:1)
if set_scaleheight
A = get(gca,'position');          % gca points at the second one
A(1,4) = A(1,4) *scaleheight;              % reduce the height by half
A(1,2) = A(1,2) +2*A(1,4)/scaleheight*(1-scaleheight);         % change the vertical position
set(gca,'position',A);     
end

end