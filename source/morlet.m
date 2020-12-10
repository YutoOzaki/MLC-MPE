classdef morlet < handle
    properties
    end
    
    methods
        function obj = morlet()
        end
    end
    
    methods(Static)
        function w = wavefun(t, w0)
            w = pi^(-0.25) .* exp(1i.*w0.*t) .* exp(-t.^2./2);
        end
        
        function w = scaledwave(s, dt, n_i, n, w0)
            % n: translation from origin
            
            w = (dt/s).^(0.5) .* morlet.wavefun((n_i - n).*dt./s, w0);
        end
        
        function g = spectrumfun(omega, w0, s)
            H = omega > 0;
            
            A = bsxfun(@times, s, omega);
            
            g = pi.^(-0.25) .* bsxfun(@times, H, exp(-(A - w0).^2./2));
        end
        
        function g = scaledspectrum(omega, w0, s, dt)
            g = bsxfun(@times, (2*pi.*s/dt).^(0.5), morlet.spectrumfun(omega, w0, s));
        end
        
        function W = transform(x, g)
            W = ifft(bsxfun(@times, fft(x), conj(g)));
        end
    end
end

%{
dt = 0.01;
t = -4:dt:4;

w0 = 6;

w = morlet.wavefun(t, w0);

mean = sum(w)/numel(w);

figure(1); 
plot(t, real(w)); hold on;
plot(t, imag(w), 'm-.'); hold off;
title(sprintf('Admissibility: mean(re), mean(im) = %e, %e\n', real(mean), imag(mean)));
%}

%{
K = 6;
S = [0.5 1 2 2 2 1];

for k=1:K
    s = S(k);

    dt = 0.01;
    t = -20:dt:20;

    n_i = 0:(numel(t) - 1);
    N = numel(n_i);

    w0 = 6;

    n = floor(N/(K + 1)) * k;
    w = morlet.scaledwave(s, dt, n_i, n, w0);

    mean = sum(w)/numel(w);

    figure(1); 
    subplot(K, 1, k);
    plot((t - t(n))./s, real(w)); hold on;
    plot((t - t(n))./s, imag(w), 'm-.'); hold off;
    title(sprintf('Admissibility: re(E[w]), im(E[w]), n, t, s = %e, %e, %d, %2.2f, %1.2f\n', real(mean), imag(mean), n, t(n), s));
end
%}

%{
dt = 0.1;
t = -20:dt:20;

n_i = 0:(numel(t) - 1);
N = numel(n_i);

w0 = 6;
s = 1;

n = ceil(N/2);
w1 = morlet.scaledwave(s, dt, n_i, n, w0);
w2 = morlet.wavefun(t, w0);

r1 = sum(real(w1) - real(w2)).^2;
r2 = sum(imag(w1) - imag(w2)).^2;
%}

%{
dt = 0.01;
t = -4:dt:4;

n_i = 0:(numel(t) - 1);
N = numel(n_i);

w0 = 6;
s = 2*rand;

k = [0:floor(N/2) -((floor(N/2) + 1):(N - 1))];
omega = 2.*pi.*k./(N*dt);
g = morlet.scaledspectrum(omega, w0, s, dt);

energy = sum(abs(g).^2)/N;

figure(1); 
subplot(2,1,1); plot(omega, g);
title(sprintf('Admissibility: energy = %1.3f, s = %1.2f\n', energy, s)); xlim([0 200]);
subplot(2,1,2); plot(g);
%}

%{
dt = 0.01;
N = 1023 + randi(2);

n_i = 0:(N - 1);

w0 = 6;
s = 2*rand;

k = [0:floor(N/2) -((floor(N/2) + 1):(N - 1))];
omega = 2.*pi.*k./(N*dt);
g = morlet.scaledspectrum(omega, w0, s, dt);

energy = sum(abs(g).^2)/N;

figure(1); 
subplot(2,1,1); plot(omega, g);
title(sprintf('Admissibility: energy = %1.3f, s = %1.2f\n', energy, s)); xlim([0 200]);
subplot(2,1,2); plot(g);
%}

%{
N = 1023 + randi(2);

dt = 0.01;
sigma = 2 * rand;
J = 3000;

w0 = 6;
k = [0:floor(N/2) -((floor(N/2) + 1):(N - 1))];
omega = 2.*pi.*k./(N*dt);

P_W = zeros(J, N);
P_hist = zeros(J, 1);
P_X1 = zeros(J, 1);
P_X2 = zeros(J, 1);

for j=1:J
    x = rand(1, N);
    x = x - mean(x);
    x = sqrt(sigma) .* x./sqrt(var(x));

    X = abs(fft(x)./N).^2;
    P_X1(j) = sum(X);

    if mod(N, 2) == 0
        X(2:N/2) = 2.*X(2:N/2);
        X(N/2+2:end) = 0;
    else
        X(2:ceil(N/2)) = 2.*X(2:ceil(N/2));
        X(ceil(N/2)+1:end) = 0;
    end
    
    P_X2(j) = sum(X);

    s = rand;
    g = morlet.scaledspectrum(omega, w0, s, dt);
    P_g = sum(abs(conj(g)).^2)/N;

    W = morlet.transform(x, g);
    P_W(j, :) = abs(W).^2;
    
    P_hist(j) = sum(sum(P_W(1:j, :)))/(j*N);

    if J < 7
        figure(1); subplot(J, 1, j); plot(conj(g)); title(sprintf('Expected energy %1.5f', P_g));
    end
end

figure(2);
subplot(2,1,1);
plot(P_hist); hold on; 
plot(sigma.*ones(J, 1), '-.m'); hold on;
plot(P_X1, '--b'); hold on;
plot(P_X2, '--g'); hold off;
subplot(2,1,2);
plot(P_hist(end-99:end)); hold on; 
plot(sigma.*ones(100, 1), '-.m'); hold on;
plot(P_X1(end-99:end), '--b'); hold on;
plot(P_X2(end-99:end), '--g'); hold off;
%}