function metric_value_cons = fi_fun(x,md)
    x_pure = x';
    
    x_0_oc = [1;1];
    time.T = 6;
    time.t = 1;

    mat_perf.Q = 0.5*eye(2); 
    mat_perf.R = 2 * eye(2);
    
    [~,x_rec] = sys_real_u(x_0_oc,time,md,x_pure,mat_perf);
    metric_value_cons = x_rec.^2-0.7^2;
    metric_value_cons = metric_value_cons';
end