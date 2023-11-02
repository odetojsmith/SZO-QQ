obj_value_es = out.simout.Data;
obj_value_es = squeeze(obj_value_es);
time_es = out.simout.Time*0.016;
save('obj_value_es.mat','obj_value_es');
save('time_es.mat','time_es');


figure(1); 
plot(time,state_value(1,:),'.-','Color','y','LineWidth',1); 
hold on;
plot(time,state_value(2,:),'.-','Color','b','LineWidth',1); 
set(gca,'Fontname','latex','linewidth',1.5,'FontSize',25,'TickLabelInterpreter','latex'); 
xlabel('$t$[s]','Interpreter','latex','FontSize',25);
ylabel('\textrm{Value of the decision variable}','Interpreter','latex','FontSize',25);

legend({'Algorithm 1','Log-barrier method'},'Interpreter','latex','Location','northeast');
