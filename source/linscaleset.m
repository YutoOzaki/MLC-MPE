function [s, f] = linscaleset(N, fs, w0)
    df = fs/N;
    
    LP = df;
    HP = fs/2;
    
    f = LP:df:HP;
    T = 1./f;
    
    s = T.*(w0 + sqrt(2 + w0^2))./(4*pi);
end