function fi = fi_fun(x)
    f1 = 0.5-(x(1)+0.5)^2-(x(2)-0.5)^2;
    f2 = x(1)^2-x(2);
    f3 = x(2)-1;

    fi = [f1,f2,f3];
end