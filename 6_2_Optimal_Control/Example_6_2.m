clear all
%% Initialization


%%%%%%%%%%%%%% Computation of the initial point
mat_perf.Q = 5*eye(2); 
mat_perf.R = 2 * eye(2);
md.A = [1.1 1;-0.5 1.1];
md.B = eye(2);
md_p.A = [1.1 1;-0.5 1.1];
md_p.B = eye(2);
time.T = 6;
time.t = 1;
up_bound = 1.2;
x_0_oc = [1;1];
[K_0,~,~] = dlqr(md_p.A,md_p.B,mat_perf.Q,mat_perf.R,0);
K_0 = -K_0;
[cost,x_rec,u_rec] = sys_real_K(x_0_oc,time,md,K_0,mat_perf);
[cost,x_rec] = sys_real_u(x_0_oc,time,md,u_rec,mat_perf);
x_0_original = u_rec;

mat_perf.Q = 0.5*eye(2); 
%%%%%%%%%%%%% Optimal point calculation (norminal)
QCell = repmat({mat_perf.Q}, 1, time.T);
RCell = repmat({mat_perf.R}, 1, time.T);
BigQ = blkdiag(QCell{:})
BigR = blkdiag(RCell{:})
u_opt = sdpvar(time.T * size(mat_perf.R ,1),1);
x_opt = sdpvar(time.T * size(mat_perf.Q ,1),1);
obj = u_opt' * BigR * u_opt + x_opt' * BigQ * x_opt;
cons = [];
for i = 1 : time.T
    x_new = x_opt((i-1) * size(mat_perf.Q ,1)+1 : i * size(mat_perf.Q ,1));
    u_temp = u_opt((i-1) * size(mat_perf.R ,1)+1 : i * size(mat_perf.R ,1));
    if i == 1
        x_old = x_0_oc;
    else
        x_old = x_opt((i-2) * size(mat_perf.Q ,1)+1 : (i-1) * size(mat_perf.Q ,1));
    end
    cons = [cons x_new == sys_dynamic_norminal(x_old,u_temp,md) x_new(1)^2<=0.7^2 x_new(2)^2<=0.7^2];
    
end
ops = sdpsettings('solver','mosek','verbose',0);
opt = optimize(cons,obj,ops);
cost_opt = value(obj);
x_opt = value(x_opt);
u_opt = value(u_opt);
[cost_ini,x_rec_ini] = sys_real_u(x_0_oc,time,md,u_opt,mat_perf);

%%%%%%%%%%%%% Optimal point calculation (norminal)
opti = casadi.Opti();
u_op = opti.variable(size(u_opt,1),1);
x_op = opti.variable(size(x_opt,1),1);
obj = u_op' * BigR * u_op + x_op' * BigQ * x_op;
for i = 1 : time.T
    x_new = x_op((i-1) * size(mat_perf.Q ,1)+1 : i * size(mat_perf.Q ,1));
    u_temp = u_op((i-1) * size(mat_perf.R ,1)+1 : i * size(mat_perf.R ,1));
    if i == 1
        x_old = x_0_oc;
    else
        x_old = x_op((i-2) * size(mat_perf.Q ,1)+1 : (i-1) * size(mat_perf.Q ,1));
    end
    opti.subject_to(x_new == sys_dynamic(x_old,u_temp,md));
    opti.subject_to(x_new(1)^2<=0.7^2);
    opti.subject_to(x_new(2)^2<=0.7^2);  
end
opti.set_initial(x_op,x_rec);
opti.set_initial(u_op,u_rec);
opti.minimize(obj);
opti.solver('ipopt');
sol = opti.solve();
sol.value(x_op)
sol.value(obj)
%%%%%%%%%%%%% 
dim = length(u_rec);        % Dimension of the problem
m = length(u_rec);          % Number of the constraints
M = 20;          % Smoothness
nu = 0.0000001;  % Fixed step length for gradient approximation
mu = 0.0001;    % Quadratic penalty coeff
xi = 0.003;  % Terminal conditions

%% Compute the optimizer
x_0 = [f_0(x_0_original)+0.1;x_0_original];       
x = x_0;
x_old = x_0+1;                  
x_diff = x-x_old;
i = 1;
while max(abs(x_diff(:))) >= xi
    grad_cons = [];
    CON = [];
    x(1) = x(1)+ 1/i;
    x_old = x;
    x_progress = sdpvar(dim+1,1);
    [metric_value_cons] = func_realization(x,md);
    for j = 1 : dim+1
        x_copy = x;
        x_copy(j) = x_copy(j)+nu;
        [metric_value_cons_pert]= func_realization(x_copy,md);
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
x_0 = [f_0(x_0_original)+0.1;x_0_original];       
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
    [metric_value_cons] = func_realization(x,md);
    for j = 1 : dim+1
        x_copy = x;
        x_copy(j) = x_copy(j)+nu;
        [metric_value_cons_pert]= func_realization(x_copy,md);
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
    
    f_0(x(2:end))
    cost = [cost f_0(x(2:end))];
    x_diff = x-x_old;
    i = i+1;
end
toc



function y=f_0(x)
%%%%%%% Here enter the evaluation of the objective function
    x_0_oc = [1;1];
    time.T = 6;
    time.t = 1;
    md.A = [1.1 1;-0.5 1.1];
    md.B = eye(2);
    mat_perf.Q = 0.5*eye(2); 
    mat_perf.R = 2 * eye(2);
    [y,~] = sys_real_u(x_0_oc,time,md,x,mat_perf);
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



function [metric_value_cons] = func_realization(x,md)
    x_pure = x(2:end);
    
    x_0_oc = [1;1];
    time.T = 6;
    time.t = 1;

    mat_perf.Q = 0.5*eye(2); 
    mat_perf.R = 2 * eye(2);
    
    [cost,x_rec] = sys_real_u(x_0_oc,time,md,x_pure,mat_perf);
    metric_value_cons = [cost-x(1); x_rec.^2-0.7^2];
end


