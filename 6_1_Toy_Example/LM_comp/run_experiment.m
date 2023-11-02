function [cost,timing,X_diff] = run_experiment(L,M)

dim = 2;        % Dimension of the problem
m = 3;          % Number of the constraints
nu = 0.000001;  % Fixed step length for gradient approximation
mu = 0.0001;    % Quadratic penalty coeff
x_0 = [0.9;0.9];               % Initial point
N = 800;

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
X_diff = [];
while i <= N
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

while i <= N
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
    x_rec = [x_rec x];

    timing = [timing toc];
    cost = [cost f_0(x(2:end))];
    x_diff = x-x_old;
    i = i+1;
    X_diff = [X_diff max(abs(x_diff(:)))];
end
x = x_rec;
X_diff = [X_diff 0];











