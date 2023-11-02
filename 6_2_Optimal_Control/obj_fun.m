function y = obj_fun(x)
    x = x';
    x_0_oc = [1;1];
    time.T = 6;
    time.t = 1;
    md.A = [1.1 1;-0.5 1.1];
    md.B = eye(2);
    mat_perf.Q = 0.5*eye(2); 
    mat_perf.R = 2 * eye(2);
    [y,~] = sys_real_u(x_0_oc,time,md,x,mat_perf);
end