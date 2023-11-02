clear all;
close all;
load('COST_L.mat')
load('COST_L_log.mat')
load('COST_M.mat')
load('COST_M_log.mat')

load('TIMING_L.mat')
load('TIMING_L_log.mat')
load('TIMING_M.mat')
load('TIMING_M_log.mat')

load('X_DIFF_L.mat')
load('X_DIFF_M.mat')

load("cost.mat")

load("cost_record_iter_safeopt.mat")
load("n_sample_log.mat")
load("n_samples.mat")
load("cost_log.mat")

%% plot of the cost (M)
figure(1)
hold on;
scale_log = 1 : 2000 : 600000;
colorset = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30'};
TIMING_M_log = TIMING_M_log(:,scale_log);
COST_M_log = COST_M_log(:,scale_log);
for i  = 1 : size(COST_L,1)
    plot(TIMING_M(i,:),COST_M(i,:),'Color',colorset{i},'LineWidth',4);
end

for i  = 1 : size(COST_L,1)
    semilogy(TIMING_M_log(i,:),COST_M_log(i,:),':','Color',colorset{i},'LineWidth',4); 
end

xlim([0 25])
ylim([10^(-3) 1])
xlabel('$\mathrm{computation}$ $\mathrm{time}$  [s]','Interpreter','latex','FontSize',25);
ylabel('$f_0(x)$','Interpreter','latex','FontSize',25);
set(gca,'XTick', 0:5:25)
set(gca, 'YScale', 'log')
set(gca,'Fontname','latex','linewidth',1.5,'FontSize',25,'TickLabelInterpreter','latex'); 
legend({'SZO-QQ: $M=3, L=5$','SZO-QQ: $M=6, L_0=5$','SZO-QQ: $M=10, L=5$','SZO-QQ: $M=15, L=5$','SZO-QQ: $M=20, L=5$','0-LBM: $M=3, L=5$','0-LBM: $M=6, L=5$','0-LBM: $M=10, L=5$','0-LBM: $M=15, L=5$','0-LBM: $M=20, L=5$'},'FontSize',20,'Interpreter','latex','Location','northeast');

figure(2)
hold on;
TIMING_L_log = TIMING_L_log(:,scale_log);
COST_L_log = COST_L_log(:,scale_log);
for i  = 1 : size(COST_L,1)
    plot(TIMING_L(i,:),COST_L(i,:),'Color',colorset{i},'LineWidth',4);
end

for i  = 1 : size(COST_L,1)
    semilogy(TIMING_L_log(i,:),COST_L_log(i,:),':','Color',colorset{i},'LineWidth',4); 
end
xlabel('$\mathrm{computation}$ $\mathrm{time}$  [s]','Interpreter','latex','FontSize',25);
ylabel('$f_0(x)$','Interpreter','latex','FontSize',25);
set(gca,'XTick', 0:5:25)
set(gca, 'YScale', 'log')
set(gca,'Fontname','latex','linewidth',1.5,'FontSize',25,'TickLabelInterpreter','latex'); 
legend({'SZO-QQ: $M=3, L=5$','SZO-QQ: $M=3, L=8$','SZO-QQ: $M=3, L=10$','SZO-QQ: $M=3, L=20$','SZO-QQ: $M=3, L=40$','0-LBM: $M=3, L=5$','0-LBM: $M=3, L=8$','0-LBM: $M=3, L=10$','0-LBM: $M=3, L=20$','0-LBM: $M=3, L=40$'},'FontSize',20,'Interpreter','latex','Location','northeast');
xlim([0 18])
ylim([10^(-3) 1])


% figure(3)
% hold on;
% TIMING_L_log = TIMING_L_log(:,scale_log);
% COST_L_log = COST_L_log(:,scale_log);
% for i  = 1 : size(COST_L,1)
%     plot(TIMING_L(i,:),COST_L(i,:),'Color',colorset{i},'LineWidth',4);
% end
% 
% for i  = 1 : size(COST_L,1)
%     semilogy(TIMING_L_log(i,:),COST_L_log(i,:),':','Color',colorset{i},'LineWidth',4); 
% end
% xlabel('$\mathrm{computation}$ $\mathrm{time}$  [s]','Interpreter','latex','FontSize',25);
% ylabel('$f_0(x)$','Interpreter','latex','FontSize',25);
% set(gca,'XTick', 0:5:25)
% set(gca, 'YScale', 'log')
% set(gca,'Fontname','latex','linewidth',1.5,'FontSize',25,'TickLabelInterpreter','latex'); 
% legend({'SZO-QQ: $M=3, L=5$','SZO-QQ: $M=3, L=8$','SZO-QQ: $M=3, L=10$','SZO-QQ: $M=3, L=20$','SZO-QQ: $M=3, L=40$','0-LBM: $M=3, L=5$','0-LBM: $M=3, L=8$','0-LBM: $M=3, L=10$','0-LBM: $M=3, L=20$','0-LBM: $M=3, L=40$'},'FontSize',20,'Interpreter','latex','Location','northeast');
% xlim([0 18])
% ylim([10^(-3) 1])

% load('obj_value_es.mat')
% load('cost_record_iter_safeopt.mat')
% figure(4); 
% n_natural = 1:length(obj_value_es);
% n_natural_2 = 1:length(cost_record_iter_safeopt);
% semilogy(n_samples,cost,'.-','Color',[0.4940 0.1840 0.5560],'LineWidth',4); 
% hold on;
% semilogy(n_sample_log,cost_log,'.-','Color','r','LineWidth',4); 
% semilogy(n_natural,obj_value_es','.-','Color',[0.9290 0.6940 0.1250],'LineWidth',4); 
% semilogy(n_natural_2,cost_record_iter_safeopt,'.-','Color',[0.4660 0.6740 0.1880],'LineWidth',4)
% 
% 
% 
% xlabel('$\mathrm{computation}$ $\mathrm{time}$  [s]','Interpreter','latex','FontSize',25);
% ylabel('$f_0(x)$','Interpreter','latex','FontSize',25);
% xlim([0 6])
% ylim([1*10^(-7) 1.5])
% set(gca,'XTick', 0:1:6)
% set(gca,'Fontname','latex','linewidth',1.5,'FontSize',25,'TickLabelInterpreter','latex'); 
% % set(gca,'TickLength',[0.025, 0.01])
% %plot(timing(end),cost(end),'.','Color','g','MarkerSize',50)
% 
% legend({'SZO-FW','SZO-QQ','s0-LBM','Extremum Seeking','SafeOptSwarm'},'Interpreter','latex','Location','northeast');
