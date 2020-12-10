function [W, z, z_q] = multilayeredceps(F, f_k, q_k, V, gam_0, gam, g)
    fs = F(2) + F(end);
    N = size(V, 1);
    
    %% setup filters
    [W_f, W_q] = filteringmat(F, N, f_k, q_k);
    
    %% 0th-layer
    z = g(abs(V./N), gam_0);
    
    %% subsequent layers
    for i=1:length(gam)
        if mod(i, 2) == 1
            h = W_q*fft(z, N, 1);
            I = h <= 0;
            h(I) = 0;
            z_q = g(h, gam(i));
        else
            h = W_f*fft(z_q, N, 1)./N;
            I = h <= 0;
            h(I) = 0;
            z = g(h, gam(i));
        end
    end
    
    %% de-shape
    T_q = (0:(N-1))./fs;
    n_q = round((1./T_q).*(N/fs));
    n_q(1) = 1;
    
    W = z.*z_q(n_q, :);
end

function [W_f, W_q] = filteringmat(F, N, f_k, q_k)
    fs = F(2) + F(end);
    
    [~, I] = sort(abs(F - f_k), 'asc');
    k_c = I(1);
    
    dt = 1/fs;
    q_c = floor((1/q_k)/dt);
    
    W_f = ones(N, 1);
    W_f(1:(k_c-1)) = 0;
    W_f((N-k_c):N) = 0;
    W_f = diag(W_f);
    
    W_q = ones(N, 1);
    W_q(1:(q_c-1)) = 0;
    W_q((N-q_c):N) = 0;
    W_q = diag(W_q);
end