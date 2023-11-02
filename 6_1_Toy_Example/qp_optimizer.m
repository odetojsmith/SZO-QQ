clear all
dim = 2;        % Dimension of the problem
m = 3;          % Number of the constraints
M = 0.8;          % Smoothness
nu = 0.000001;  % Fixed step length for gradient approximation
mu = 0.0001;    % Quadratic penalty coeff
xi = 0.000001;  % Terminal conditions
x_0 = [0.9;0.9];               % Initial point
L = 0.8;


%% Compute the optimizer
x_0 = [f_0(x_0)+0.1;x_0]; 
x_rec = x_0; 
x = x_0;
x_old = x_0+1;                  
x_diff = x-x_old;
i = 1;
alpha = sqrt(dim) * M/2;
Lambda = 3;
C = 4 * Lambda * (alpha*m+m*L+4*m*M);
ETA = [1 0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001 0.0005];
XI = ETA/C;
pointer = 1;
time_xi = [];
while max(abs(x_diff(:))) >= XI(end)
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
x_0 = [0.9;0.9];
x_0 = [f_0(x_0)+0.1;x_0];
x = x_0;
x_old = x_0+1;
x_diff = x-x_old;
i = 1;
timing = 0;
cost = f_0(x(2:end));
count = 0;
n_samples = 1;
tic

while max(abs(x_diff(:))) >= XI(end)
    grad_cons = [];
    count = count+1;

    CON = [];
    x_old = x;
    x_progress = sdpvar(dim+1,1);
    t_sample = 0;
    t_1 = toc;
    [metric_value_cons] = func_realization(x,m);
    t_2 = toc;
    t_sample = t_sample + t_2-t_1;
    for j = 1 : dim+1
        x_copy = x;
        x_copy(j) = x_copy(j)+nu;
        t_1 = toc;
        [metric_value_cons_pert]= func_realization(x_copy,m);
        t_2 = toc;
        t_sample = t_sample + t_2-t_1;
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

    x = P(x,grad_cons,metric_value_cons);
    x_rec = [x_rec x];

    timing = [timing toc-t_sample];
    cost = [cost f_0(x(2:end))];
    x_diff = x-x_old;
    i = i+1;
    if max(abs(x_diff(:))) < XI(pointer)
        time_xi = [time_xi toc];
        pointer= pointer+1;
    end
    n_samples = [n_samples 1+3*count];
end
x = x_rec;
save('x.mat','x');
save('cost.mat','cost');
save('timing.mat','timing');
save('n_samples.mat','n_samples');

save('ETA.mat','ETA');
save('time_xi.mat','time_xi');







function y=f_0(x)
%%%%%%% Here enter the evaluation of the objective function
    y = x(2)+0.1*x(1)^2;
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
    if k == 1
        y = 0.5-(x(1)+0.5)^2-(x(2)-0.5)^2;
    end
    if k == 2
        y = x(2)-1;
    end
    if k == 3
        y = x(1)^2-x(2);
    end
end

function [metric_value_cons] = func_realization(x,m)
    x_pure = x(2:end);
    metric_value_cons = f_0(x_pure)-x(1);
    for i = 1 : m
        metric_value_cons = [metric_value_cons;f_c(x_pure,i)];
    end
end


