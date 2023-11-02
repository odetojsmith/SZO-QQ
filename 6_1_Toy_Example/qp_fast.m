clear all
dim = 2;
m = 4;
M = 3;
nu = 0.000001;
mu = 0.0001;
xi = 0.0005;

x_0 = [1.99;1.6];
x_0 = [f_0(x_0)+0.1;x_0];
x = x_0;
x_old = x_0+1;
x_diff = x-x_old;
i = 1;
tic
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
    obj = [1 zeros(1,dim)] * x_progress + mu * (x_progress-x(:))' * (x_progress-x(:));
    
    for j = 1 : m+1
        CON_j = metric_value_cons(j) + grad_cons(:,j)' * ((x_progress-x)) + 2*M * (x_progress-x)'* (x_progress-x)<=0;
        CON = [CON, CON_j];
    end

    x_out = opt_fast(metric_value_cons,x,grad_cons,M,dim,m+1);
    x_diff = x-x_old;
    i = i+1
end
toc

cost  = f_0(x(2:end))

function y=f_0(x)
y = norm(x + [-1;1]).^2;
end

function y=f_c(x,k)
    if k == 1
        y = x(1)-2;
    end
    if k == 2
        y = -x(1);
    end
    if k == 3
        y = x(2)-2;
    end
    if k == 4
        y = -x(2);
    end
end

function [metric_value_cons] = func_realization(x,m)
    x_pure = x(2:end);
    metric_value_cons = f_0(x_pure)-x(1);
    for i = 1 : m
        metric_value_cons = [metric_value_cons;f_c(x_pure,i)];
    end
end


