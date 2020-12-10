function c = pitchtracking(x, W_q, eps_t, eps_c, col_I)
    c = zeros(1, size(x, 2));
    T = unique(x(1, :));
    
    K = 0;
    UB = normpdf(0);
    
    W_t = W_q;
    W_t(W_t ~= 0) = 1;
    N_t = sum(W_t);
    N_t(N_t == 0) = [];
    
    idx_s = 1;
    
    for t=1:(length(T) - 1)
        %% find data points at the current time index
        idx_e = idx_s + N_t(t) - 1;
        p_idx = idx_s:idx_e;
        idx_s = idx_e + 1;

        %% prioritize the beginning of tracking with saliency
        E = W_q(W_q(:, col_I(t)) ~= 0, col_I(t));
        [~, I] = sort(E, 'desc');
        
        %% perform tracking per note
        for j=1:length(I)
            t_step = 1;
            idx_et = idx_e;
            continuetracking = true;
            
           %% assign cluster if not yet
            i = p_idx(I(j));
            if c(i) == 0
                K = K + 1;
                c(i) = K;
            end
                
            while continuetracking
                t_next = T(t + t_step);
                
                if (t_next - T(t)) > eps_t
                    break;
                else
                  %% find data points at the next time index
                    idx_st = idx_et + 1;
                    idx_et = idx_st + N_t(t + t_step) - 1;
                    q_idx = idx_st:idx_et;
                    
                  %% candidate's cents
                    c_q = x(2, q_idx);
                    
                  %% sort candidate's with saliency
                    E_t = W_q(W_q(:, col_I(t + t_step)) ~= 0, col_I(t + t_step));
                    [~, I_t] = sort(E_t, 'desc');
                    
                  %% pick up candidate notes
                    r = normpdf(abs(x(2, i) - c_q)./eps_c);
                    p = rand*UB;
                    
                    for l=1:length(I_t)
                        if p < r(I_t(l))
                            if c(q_idx(I_t(l))) == 0
                                c(q_idx(I_t(l))) = c(i);

                                continuetracking = false;
                            end
                        end
                    end
                    
                  %% update time step
                    if (t + t_step) ~= length(T)
                        t_step = t_step + 1;
                    else
                        continuetracking = false;
                    end
                end
            end
        end
    end
end