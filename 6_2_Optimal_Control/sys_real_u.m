function [cost,x_rec] = sys_real_u(x_0,time,md,u_vec,mat_perf)
    cost = 0;
    x = x_0;
    x_rec = [];
    dim_u = size(mat_perf.Q,1);
    for i = 1:time.T
       u = u_vec((dim_u) * (i-1) + 1 : dim_u * i);
       x = sys_dynamic(x,u,md); 
       cost = cost + u' * mat_perf.R * u + x' * mat_perf.Q * x;
       x_rec = [x_rec; x];
    end
end