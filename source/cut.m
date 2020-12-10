function x_t = cut(x, fs, x_start, x_end)
    dt = 1/fs;
    
    if (dt*(length(x)-1)) < x_end
        x_end = dt*(length(x)-1);
    end
    
    t_s = floor(x_start/dt) + 1;
    t_e = ceil(x_end/dt) + 1;
    
    x_t = x(t_s:t_e);
end