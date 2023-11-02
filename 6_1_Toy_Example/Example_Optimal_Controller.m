clear all
%% Initialization
dim = 4;        % Dimension of the problem
m = 4;          % Number of the constraints
M = 3;          % Smoothness
nu = 0.000001;  % Fixed step length for gradient approximation
mu = 0.0001;    % Quadratic penalty coeff
xi = 0.000001;  % Terminal conditions

%%%%%%%%%%%%%% Computation of the initial point
mat_perf.Q = 0.1*eye(2);
mat_perf.R = eye(2);
md.A = [1.1 1;-0.5 1.1];
md.B = eye(2);
md_p.A = [1.5 1;-0.5 1.5];
md_p.B = eye(2);
time.T = 15;
time.t = 1;
up_bound = 1.2;
x_0_oc = [1;1];
[K_0,~,~] = dlqr(md_p.A,md_p.B,mat_perf.Q,mat_perf.R,0);
[K_opt,~,~] = dlqr(md.A,md.B,mat_perf.Q,mat_perf.R,0);
K_0 = -K_0;
[cost,x_quest] = sys_realization(x_0_oc,time,md,K_0,mat_perf);
[cost_opt,x_quest_opt] = sys_realization(x_0_oc,time,md,-K_opt,mat_perf);
x_0 = K_0(:);               % Initial point

%% Compute the optimizer
x_0 = [f_0(x_0)+0.1;x_0];       
x = x_0;
x_old = x_0+1;                  
x_diff = x-x_old;
i = 1;
while max(abs(x_diff(:))) >= xi
    grad_cons = [];
    CON = [];
    x_old = x;
    x_progress = sdpvar(dim+1,1);
    [metric_value_cons] = func_realization(x,m);
    for j = 1 : dim+1
        x_copy = x;
        x_copy(j) = x_copy(j)+nu;
        [metric_value_cons_pert]= func_realization(x_copy,m);
        grad_cons_j = (metric_value_cons_pert-metric_value_cons)/nu;
        grad_cons = [grad_cons;grad_cons_j'];
    end
    if i == 1
        x_var = sdpvar(dim+1,1);
        grad_var = sdpvar(dim+1,m+1);
        metric_var = sdpvar(m+1,1);
        obj = [1 zeros(1,dim)] * x_progress + mu * (x_progress-x_var)' * (x_progress-x_var);

        for j = 1 : m+1
            CON_j = metric_var(j) + grad_var(:,j)' * ((x_progress-x_var)) + 2*M * (x_progress-x_var)'* (x_progress-x_var)<=0;
            CON = [CON, CON_j];
        end

        P = optimizer(CON,obj,[],{x_var,grad_var,metric_var},x_progress);
        x = P(x,grad_cons,metric_value_cons);
        break;
    end
    x = P(x,grad_cons,metric_value_cons);
    
    x_diff = x-x_old;
    i = i+1;
end

%% Start the iteration
x_0 = K_0(:);
x_0 = [f_0(x_0)+0.1;x_0];
x = x_0;
x_old = x_0+1;
x_diff = x-x_old;
i = 1;
timing = 0;
cost = f_0(x(2:end));

while max(abs(x_diff(:))) >= xi
    grad_cons = [];

    CON = [];
    x_old = x;
    x_progress = sdpvar(dim+1,1);
    [metric_value_cons] = func_realization(x,m);
    for j = 1 : dim+1
        x_copy = x;
        x_copy(j) = x_copy(j)+nu;
        [metric_value_cons_pert]= func_realization(x_copy,m);
        grad_cons_j = (metric_value_cons_pert-metric_value_cons)/nu;
        grad_cons = [grad_cons;grad_cons_j'];
    end
    
    x_var = sdpvar(dim+1,1);
    grad_var = sdpvar(dim+1,m+1);
    metric_var = sdpvar(m+1,1);
    obj = [1 zeros(1,dim)] * x_progress + mu * (x_progress-x_var)' * (x_progress-x_var);

    for j = 1 : m+1
        CON_j = metric_var(j) + grad_var(:,j)' * ((x_progress-x_var)) + 2*M * (x_progress-x_var)'* (x_progress-x_var)<=0;
        CON = [CON, CON_j];
    end
    if i == 1
        tic
    end
    x = P(x,grad_cons,metric_value_cons);
    
    timing = [timing toc];
    cost = [cost f_0(x(2:end))];
    x_diff = x-x_old;
    i = i+1;
end



function y=f_0(x)
%%%%%%% Here enter the evaluation of the objective function
    x_0_oc = [1;1];
    time.T = 15;
    time.t = 1;
    md.A = [1.1 1;-0.5 1.1];
    md.B = eye(2);
    mat_perf.Q = 0.1*eye(2);
    mat_perf.R = eye(2);
    
    x_oc_dim = sqrt(length(x));
    K = reshape(x,[x_oc_dim,x_oc_dim]);
    [y,~] = sys_realization(x_0_oc,time,md,K,mat_perf);
end

% function y=f_c(x,k)
%     if k == 1
%         y = x(1)-2;
%     end
%     if k == 2
%         y = -x(1);
%     end
%     if k == 3
%         y = x(2)-2;
%     end
%     if k == 4
%         y = -x(2);
%     end
% end

function y=f_c(x,k)
%%%%%%%%% Here, enter the evaluation of the constraint functions
    x_0_oc = [1;1];
    time.T = 15;
    time.t = 1;
    md.A = [1.1 1;-0.5 1.1];
    md.B = eye(2);
    mat_perf.Q = 0.1*eye(2);
    mat_perf.R = eye(2);
    
    x_oc_dim = sqrt(length(x));
    K = reshape(x,[x_oc_dim,x_oc_dim]);
    [~,x_value] = sys_realization(x_0_oc,time,md,K,mat_perf);
    if k == 1
        y = x_value(1)-1.2;
    end
    if k == 2
        y = -x_value(1)-1.2;
    end
    if k == 3
        y = x_value(2)-1.2;
    end
    if k == 4
        y = -x_value(2)-1.2;
    end
end

function [metric_value_cons] = func_realization(x,m)
    x_pure = x(2:end);
    metric_value_cons = f_0(x_pure)-x(1);
    for i = 1 : m
        metric_value_cons = [metric_value_cons;f_c(x_pure,i)];
    end
end


