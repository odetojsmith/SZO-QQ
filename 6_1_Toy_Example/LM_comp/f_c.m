function y=f_c(x,k)
%%%%%%%%% Here, enter the evaluation of the constraint functions
    if k == 1
        y = -x(1);
    end
    if k == 2
        y = x(2)-1;
    end
    if k == 3
        y = x(1)^2-x(2);
    end
end