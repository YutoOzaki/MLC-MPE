function A4_freq = refpitch(x_t, fs, stdfreq, numsemi)
    K = 200;
    numepoch = 2000;
    n_min = 15;
    n_max = 110;
    lr = 1e-2;
    
    P = abs(fft(x_t)).^2;
    F = (0:(length(x_t)-1)) .* (fs/length(x_t));
    
    [~, I] = sort(P(1:ceil(length(x_t)/2)), 'desc');
    
    y = F(I(1:K))';
    p = P(I(1:K));
    p = p./sum(p);
    n = zeros(K, 1);
    
    A4_hist = zeros(1, numepoch);
    A4_freq = stdfreq;
    for epoch=1:numepoch
        A4_hist(epoch) = A4_freq;
        
        f = A4_freq.*2.^(((n_min:n_max) - 69)./numsemi);

        for i=1:K
            [~, I] = min((f - y(i)).^2);
            n(i) = I + n_min - 1;
        end

        yhat = A4_freq.*2.^((n - 69)./numsemi);

        d = 2.^((n - 69)./numsemi) .* p.*(yhat - y);

        A4_freq = A4_freq - lr*sum(d);
    end
end