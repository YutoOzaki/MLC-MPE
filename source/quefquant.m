function [q_q, F_q] = quefquant(q, fs, F_min, F_max, reffreq, numsemi)
    p_min = 69 + ceil(numsemi*log2(F_min/reffreq));
    p_max = 69 + floor(numsemi*log2(F_max/reffreq));
    
    assert(F_min < reffreq*2^((p_min - 69)/numsemi), '');
    assert(F_max > reffreq*2^((p_max - 69)/numsemi), '');
    
    N = size(q, 1);
    Q = fs./(0:(N-1));
    p_qidx = 69 + round(numsemi.*log2(Q./reffreq));
    
    q_q = zeros(p_max - p_min + 1, size(q, 2));
    
    i_max = N;
    for i=p_min:p_max
        q_t = find(p_qidx(1:i_max) == i);
        
        if ~isempty(q_t)
            [~, I_q] = max(q(q_t, :), [], 1);
            
            for t=1:size(q, 2)
                q_q(i - p_min + 1, t) = q(q_t(I_q(t)), t);
            end
            
            i_max = q_t(1) - 1;
        end
    end
    
    F_q = reffreq.*2.^(((p_min:p_max) - 69)./numsemi);
end