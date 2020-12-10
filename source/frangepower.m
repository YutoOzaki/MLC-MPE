function [f_lb, f_ub, S] = frangepower(x_t, fs, w0, dt)
    A4 = 440;
    dc = 10;
    beta = 0.1;
    
    eps_u = 0.7;
    eps_l = 0.075;
    
    f_min = 25;
    f_max = 4000;
    
    [s, f] = logscaleset(A4, dc, w0);
    
    omega = angfreqset(x_t, fs);
    
    gm = morlet.scaledspectrum(omega, w0, s, 1/fs);
    
    V = morlet.transform(x_t, gm);
    
    dn = round(fs*dt);
    V = V(1:dn:end, :);
    
    E = max((abs(V).^2).^beta, [], 1);
    E = E./sum(E);
    E = cumsum(E);
    
    I = min(find(E > eps_u));
    f_ub = min(f(I + 1), f_max);
    
    I = min(find((E - eps_l) > 0));
    f_lb = max(f(I), f_min);

    %%{  
    P = sum(abs(V).^2, 2);
    E = mean(P);
    S = P > (0.01*E);
    %}
    
    %{
    P = log(sum(abs(V).^2, 2));
    S = P > 0;
    %}
    
    %{
    P = bsxfun(@times, abs(V).^2, S);
    
    figure;
    surf(1:size(V, 1), f, log(P'), 'edgecolor', 'none');
        axis tight;
        view(0, 90);
        set(gca, 'YScale', 'log');
        grid off;
    %}
end