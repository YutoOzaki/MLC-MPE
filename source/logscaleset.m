function [s, f] = logscaleset(A4, dc, w0)
    LP = A4/(2^4);
    LC = 1200*log2(LP/A4);
    
    HP = A4*(2^4);
    HC = 1200*log2(HP/A4);
    
    f = A4 .* 2.^((LC:dc:HC)./1200);
    T = 1./f;
    
    s = T.*(w0 + sqrt(2 + w0^2))./(4*pi);
end