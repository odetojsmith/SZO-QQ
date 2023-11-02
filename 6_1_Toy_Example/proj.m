function y = proj(x,o,r)
    if norm(x-o) >= r
        y = o + (x-o)/(norm(x-o))*r;
    else
        y = x;
    end
    
end