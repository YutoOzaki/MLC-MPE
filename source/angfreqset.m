function omega = angfreqset(x_t, fs)
    N = length(x_t);
    dt = 1/fs;
    
    k = [0:floor(N/2) -((floor(N/2) + 1):(N - 1))]';
    
    omega = 2.*pi.*k./(N*dt);
end