function w = gausswindow(N, sigma)
    if nargin == 1
        sigma = 0.04;
    end
    
    t = linspace(-0.5, 0.5, N);
    
    f_w = @(x, s) 1/sqrt(2*pi*s^2) .* exp(-0.5.*(x.^2)./s^2);
    
    w = f_w(t, sigma);
end

%%
% figure; plot(fftshift(log(abs(fft(w)./N).^2)));
%