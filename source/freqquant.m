function [z_q, F_q] = freqquant(z, fs, F_min, F_max, reffreq, numsemi)
    p_min = 69 + ceil(numsemi*log2(F_min/reffreq));
    p_max = 69 + floor(numsemi*log2(F_max/reffreq));
    
    %assert(F_min < reffreq*2^((p_min - 69)/numsemi), sprintf('F_min vs. %3.3f\n', reffreq*2^((p_min - 69)/numsemi)));
    %assert(F_max > reffreq*2^((p_max - 69)/numsemi), '');
    
    N = size(z, 1);
    F = (0:(N-1)).*(fs/N);
    p_fidx = 69 + round(numsemi.*log2(F./reffreq));
    
    z_q = zeros(p_max - p_min + 1, size(z, 2));
    
    i_min = 1;
    for i=p_min:p_max
        p_t = find(p_fidx(i_min:end) == i);
        
        if ~isempty(p_t)
            p_t = p_t + i_min - 1;
            
            [~, I_f] = max(z(p_t, :), [], 1);
            
            for t=1:size(z_q, 2)
                z_q(i - p_min + 1, t) = z(p_t(I_f(t)), t);
            end
            
            i_min = p_t(end) + 1;
        end
    end
    
    F_q = reffreq.*2.^(((p_min:p_max) - 69)./numsemi);
end