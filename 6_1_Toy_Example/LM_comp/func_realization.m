function [metric_value_cons] = func_realization(x,m)
    x_pure = x(2:end);
    metric_value_cons = f_0(x_pure)-x(1);
    for i = 1 : m
        metric_value_cons = [metric_value_cons;f_c(x_pure,i)];
    end
end