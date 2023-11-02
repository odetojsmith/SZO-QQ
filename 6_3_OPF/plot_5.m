clc; clear; close all

% data in paper
% load('../data/rec_cost_LP.mat');
% load('../data/rec_max_con_LP.mat');
% load('../data/rec_timing_LP.mat');
% load('../data/rec_cost_LB.mat');
% load('../data/rec_max_con_LB.mat');
% load('../data/rec_timing_LB.mat');
% load('../data/rec_cost_QQ.mat')
% load('../data/rec_max_con_QQ.mat');
% load('../data/rec_timing_QQ.mat')
% load('../data/rec_cost_ES.mat');
% load('../data/rec_max_con_ES.mat');
% load('../data/rec_timing_ES.mat');

% data from code
load('data/SZO_LP/rec_cost_LP.mat');
load('data/SZO_LP/rec_max_con_LP.mat');
load('data/SZO_LP/rec_timing_LP.mat');
load('data/LB_SGD/rec_cost_LB.mat');
load('data/LB_SGD/rec_max_con_LB.mat');
load('data/LB_SGD/rec_timing_LB.mat');
load('data/SZO_QQ/rec_cost_QQ.mat')
load('data/SZO_QQ/rec_max_con_QQ.mat');
load('data/SZO_QQ/rec_timing_QQ.mat');
load('data/SZO_QQ/rec_cost.mat');
load('data/ES/rec_cost_ES.mat');
load('data/ES/rec_max_con_ES.mat');
load('data/ES/rec_timing_ES.mat');
load('data/cost_safeopt.mat');

cost_safeopt = -cost_safeopt;
n_sample_vec_safeopt = 1:length(cost_safeopt);




% color map
clm2 = [ 57 132 195;
         77 175 74;
         162 83 175;
         252 125 2;]/255;

%%
opti = 800.14;
rec_timing_LP(1) = 1e-2;
rec_timing_LB(1) = 1e-2;
rec_timing_QQ(1) = 1e-2;
rec_timing_ES(1) = 1e-2;


rec_timing_QQ = [rec_timing_QQ 1500];
rec_cost_QQ = [rec_cost_QQ rec_cost_QQ(end)];
rec_max_con_QQ = [rec_max_con_QQ rec_max_con_QQ(end)];


% figure(1)
% subplot(2,1,1)
% semilogx(rec_timing_LP,rec_cost_LP,'Color',clm2(1,:),'LineWidth',6); hold on
% semilogx(rec_timing_LB,rec_cost_LB,'Color',clm2(2,:),'LineWidth',6); hold on
% semilogx(rec_timing_QQ,rec_cost_QQ,'Color',clm2(3,:),'LineWidth',6); hold on
% semilogx(rec_timing_ES,rec_cost_ES,'Color',clm2(4,:),'LineWidth',6); hold on
% 
% 
% semilogx(rec_timing_QQ,ones(size(rec_timing_QQ))*opti,'--','Color','k','LineWidth',3);
% set(gca,'XTick', [1e-2 1e0 1e2 1e3]);
% set(gca,'YTick', 810:20:870)
% 
% ylim([800,882]);
% xlabel('computation time (s)','Interpreter','latex');
% ylabel('generation cost','Interpreter','latex');
% legend("SZO-LP","LB-SGD","SZO-QQ","Extremum-Seeking","Model-based",'Interpreter','latex','FontSize',26,'location','southwest');
% xlim ([0,1500]);
% set(gca,'Fontname','latex','linewidth',1.5,'FontSize',35,'TickLabelInterpreter','latex');
% grid on
% 
% ax = gca;
% ax.Position(1) = 0.13;
% ax.Position(2) = 0.63;
% ax.Position(3) = 0.85;
% ax.Position(4) = 0.36;
% 
% 
% subplot(2,1,2)
% semilogx(rec_timing_LP,rec_max_con_LP,'Color',clm2(1,:),'LineWidth',6); hold on
% semilogx(rec_timing_LB,rec_max_con_LB,'Color',clm2(2,:),'LineWidth',6); hold on
% semilogx(rec_timing_QQ,rec_max_con_QQ,'Color',clm2(3,:),'LineWidth',6); hold on
% semilogx(rec_timing_ES,rec_max_con_ES,'Color',clm2(4,:),'LineWidth',6); hold on
% 
% semilogx(rec_timing_QQ,zeros(size(rec_timing_QQ)),'--','Color','k','LineWidth',3);
% xlabel('computation time (s)','Interpreter','latex');
% ylabel('$\max_{i\geq1}f_i(x)$','Interpreter','latex');
% legend("SZO-LP","LB-SGD","SZO-QQ","Extremum-Seeking","Threshold",'Interpreter','latex','FontSize',26,'location','southwest');
% set(gca,'Fontname','latex','linewidth',1.5,'FontSize',35,'TickLabelInterpreter','latex');
% grid on
% 
% xlim ([0,1500]);
% 
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'XTick', [1e-2 1e0 1e2 1e3]);
% ax = gca;
% ax.Position(2) = 0.12;
% ax.Position(3) = 0.85;
% ax.Position(4) = 0.36;
% 
% 
% f = gcf;
% f.Position(3) = f.Position(3)*1.8;
% f.Position(4) = f.Position(4)*1.8;


%%
max_length = 3000;
rec_cost_LP_long = [rec_cost_LP rec_cost_LP(end)*ones(1,max_length-size(rec_cost_LP,2))];
rec_cost_LB_long = [rec_cost_LB rec_cost_LB(end)*ones(1,max_length-size(rec_cost_LB,2))];
rec_cost_QQ_long = [rec_cost_QQ rec_cost_QQ(end)*ones(1,max_length-size(rec_cost_QQ,2))];
rec_cost_ES_long = rec_cost_ES;

rec_con_LP_long = [rec_max_con_LP rec_max_con_LP(end)*ones(1,max_length-size(rec_max_con_LP,2))];
rec_con_LB_long = [rec_max_con_LB rec_max_con_LB(end)*ones(1,max_length-size(rec_max_con_LB,2))];
rec_con_QQ_long = [rec_max_con_QQ rec_max_con_QQ(end)*ones(1,max_length-size(rec_max_con_QQ,2))];
rec_con_ES_long = rec_max_con_ES;



figure(2)
n_sample_vec = 1 : 12 : 60000;
n_sample_vec = n_sample_vec(1:3000);
n_sample_vec_safeopt = [n_sample_vec_safeopt 1000];
cost_safeopt = [cost_safeopt; cost_safeopt(end)]';

plot(n_sample_vec(1:1777), rec_cost,'Color','b','LineWidth',6); hold on
plot(n_sample_vec, rec_cost_LB_long,'Color','r','LineWidth',6); hold on
plot(real(rec_cost_ES(1:12302)),'Color',[0.9290 0.6940 0.1250],'LineWidth',6); hold on

plot(ones(1,max(n_sample_vec)).*800.14,'--','Color','k','LineWidth',3);
set(gca,'YTick', 810:20:870)
xlim([0,13000]);
ylim([800,881]);
xlabel('number of samples','Interpreter','latex');
ylabel('generation cost','Interpreter','latex');
legend("SZO-QQ","LB-SGD","Extremum-Seeking","Reference",'Interpreter','latex','FontSize',26,'location','southeast');


set(gca,'Fontname','latex','linewidth',1.5,'FontSize',26,'TickLabelInterpreter','latex');
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'XTick', [0 2000 4000 6000 8000 10000 12000]);

grid on

ax = gca;
ax.Position(1) = 0.13;
ax.Position(3) = 0.85;
ax.Position(4) = 0.7;

f = gcf;
f.Position(3) = f.Position(3)*1.8;
f.Position(4) = f.Position(4)*0.9;



