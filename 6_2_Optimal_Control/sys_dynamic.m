function x_new = sys_dynamic(x,u,md)

    x_new = md.A * x + md.B * u +[0.1*(x(2)^2) ;0];
end