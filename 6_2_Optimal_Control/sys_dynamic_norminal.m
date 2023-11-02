function x_new = sys_dynamic_norminal(x,u,md)

    x_new = md.A * x + md.B * u;
end