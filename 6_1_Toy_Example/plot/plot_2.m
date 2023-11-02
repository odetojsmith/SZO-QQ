load("cost.mat")
load("obj_value_es.mat")
load("cost_record_iter_safeopt.mat")
load("n_sample_log.mat")
load("n_samples.mat")
load("cost_log.mat")

figure(4); 
n_natural = 1:length(obj_value_es);
n_natural_2 = 1:length(cost_record_iter_safeopt);
semilogy(n_samples,cost,'.-','Color','b','LineWidth',4); 
hold on;
semilogy(n_sample_log,cost_log,'.-','Color','r','LineWidth',4); 
semilogy(n_natural,obj_value_es','.-','Color',[0.9290 0.6940 0.1250],'LineWidth',4); 
semilogy(n_natural_2,cost_record_iter_safeopt,'.-','Color',[0.4660 0.6740 0.1880],'LineWidth',4)



xlabel('$\mathrm{number}$ $\mathrm{of}$ $\mathrm{samples}$','Interpreter','latex','FontSize',25);
ylabel('$f_0(x)$','Interpreter','latex','FontSize',25);
xlim([0 82])
ylim([1*10^(-7) 1.5])
set(gca,'XTick', 0:10:82)
set(gca,'Fontname','latex','linewidth',1.5,'FontSize',25,'TickLabelInterpreter','latex'); 
% set(gca,'TickLength',[0.025, 0.01])

legend({,'SZO-QQ','LB-SGD','Extremum Seeking','SafeOptSwarm'},'Interpreter','latex','Location','northeast');