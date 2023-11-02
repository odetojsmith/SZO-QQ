clear all
load('timing.mat')
load('timing_log.mat')
load('x.mat')
load('x_log.mat')

mu = 0.0001;
dim = 2;
m = 3;          
M = 3; 
param.mu = mu;
param.dim = dim;
param.M = M;
param.m = m;


t = [];
x_time = [];
x_time_log = [];
eta = [];
eta_log = [];
for i = 0.01: 0.01: 4.85
    t = [t i];
    vec_time = find(timing<i);
    x_time = [x_time x(:,vec_time(end))];
    vec_time_log = find(timing_log<i);
    x_time_log = [x_time_log x_log(:,vec_time_log(end))];
end
for i = 1 : length(x_time)
    if i>1 && isequal(x_time(:,i),x_time(:,i-1))
        eta = [eta eta(end)];
    else
        eta = [eta dual_calculation(param, x_time(:,i))];
    end
    if i>1 && isequal(x_time_log(:,i),x_time_log(:,i-1))
        eta_log = [eta_log eta_log(end)];
    else
        eta_log = [eta_log dual_calculation(param, x_time_log(:,i))];
    end
end
% save('t.mat','t');
% save('eta.mat','eta');
% save('eta_log.mat','eta_log');
%%
load('t.mat');
load('eta.mat');
load('eta_log.mat');
figure(1); 
semilogy(t,eta.^(-2),'.-','Color','b','LineWidth',4); 
hold on;
semilogy(t,eta_log.^(-2),'.-','Color','r','LineWidth',4); 

xlim([0 4.2])
xlabel('$\mathrm{computation}$ $\mathrm{time}$  [s]','Interpreter','latex','FontSize',25);
ylabel('$\eta^{-2}$','Interpreter','latex','FontSize',25);
set(gca,'XTick', 0:1:4.2)
set(gca,'Fontname','latex','linewidth',1.5,'FontSize',25,'TickLabelInterpreter','latex'); 
% set(gca,'TickLength',[0.025, 0.01])
legend({'Algorithm 1','Log-barrier method'},'Interpreter','latex','Location','northeast');
function y=f_0(x)
%%%%%%% Here enter the evaluation of the objective function
    y = x(2)+0.1*x(1)^2;
end


function y=f_c(x,k)
%%%%%%%%% Here, enter the evaluation of the constraint functions
    if k == 1
        y = -x(1);
    end
    if k == 2
        y = x(2)-1;
    end
    if k == 3
        y = x(1)^2-x(2);
    end
end

function grad_cons=f_c_gradient(x)
%%%%%%%%% Calculate the concrete gradient
    grad_cons = zeros(3,4);
    grad_cons(:,1) = [-1;0.02*x(2);1];
    grad_cons(:,2) = [0;-1;0];
    grad_cons(:,3) = [0;0;1];
    grad_cons(:,4) = [0;2*x(2);-1];
end

function [metric_value_cons] = func_realization(x,m)
    x_pure = x(2:end);
    metric_value_cons = f_0(x_pure)-x(1);
    for i = 1 : m
        metric_value_cons = [metric_value_cons;f_c(x_pure,i)];
    end
end

function eta = dual_calculation(param, x)
    [metric_value_cons] = func_realization(x,param.m);
    grad_cons = f_c_gradient(x);
    f_0_grad = [1;zeros(param.dim,1)];
    lambda_var = sdpvar(param.m+1,1);
    eta_var = sdpvar(1,1);
    obj = eta_var;
    CONS = [];
    CONS = [CONS norm(f_0_grad + grad_cons * lambda_var)<= eta_var];
    for i = 1:param.m+1
        CONS = [CONS -lambda_var(i) * metric_value_cons(i)<= eta_var];
    end
    CONS = [CONS lambda_var>=0];
    CONS = [CONS lambda_var<=10];
    
    ops = sdpsettings('solver','mosek','verbose',0);
    optimize(CONS,obj,ops);
    eta = value(eta_var);
end

