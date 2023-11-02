load('../cost.mat')
load('../timing.mat')

load('../cost_log.mat')
load('../timing_log.mat')
load('../obj_value_es.mat')
load('../time_es.mat')
load('../cost_record_safeopt.mat')
load('../time_record_safeopt.mat')

t = 0;
c = cost(1);
c_log = cost(1);

for i = 0.01: 0.01: 6
    t = [t i];
    vec_c = find(timing<i);
    c = [c cost(vec_c(end))];
    vec_c_log = find(timing_log<i);
    c_log = [c_log cost_log(vec_c_log(end))];
end

c = [c c(end) * ones(1,length(t)-length(c))];

figure(1); 
%semilogy(timing_fw,cost_fw,'.-','Color',[0.4940 0.1840 0.5560],'LineWidth',4); 
semilogy(t,c,'.-','Color','b','LineWidth',4); 
hold on;
semilogy(t,c_log,'.-','Color','r','LineWidth',4); 
semilogy(time_es',obj_value_es','.-','Color',[0.9290 0.6940 0.1250],'LineWidth',4); 
semilogy(time_record_safeopt,cost_record_safeopt,'.-','Color',[0.4660 0.6740 0.1880],'LineWidth',4)



xlabel('$\mathrm{computation}$ $\mathrm{time}$  [s]','Interpreter','latex','FontSize',25);
ylabel('$f_0(x)$','Interpreter','latex','FontSize',25);
xlim([0 6])
ylim([1*10^(-7) 1.5])
set(gca,'XTick', 0:1:6)
set(gca,'Fontname','latex','linewidth',1.5,'FontSize',25,'TickLabelInterpreter','latex'); 
% set(gca,'TickLength',[0.025, 0.01])
plot(timing(end),cost(end),'.','Color','g','MarkerSize',50)

legend({'SZO-QQ','s0-LBM','Extremum Seeking','SafeOptSwarm'},'Interpreter','latex','Location','northeast');
