function [cost,x_rec,u_rec] = sys_real_K(x_0,time,md,K,mat_perf)
    cost = 0;
    x = x_0;
    x_rec = [];
    u_rec = [];
    for i = 1:time.T
       u = K * x;
       x = sys_dynamic(x,u,md); 
       cost = cost + u' * mat_perf.R * u + x' * mat_perf.Q * x;
       x_rec = [x_rec;x];
       u_rec = [u_rec;u];
    end
end