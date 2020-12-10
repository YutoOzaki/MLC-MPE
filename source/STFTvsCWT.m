function STFTvsCWT
    %% parameters
    fs_d = 22050;
    
    N = 4096;
    M = round(0.98*N);
    w = gausswindow(N);
    
    w0 = 64;
    dt = 0.005;
    
    %% test cases
    dataset = mpedata();
    
    %% run
    for n=1:length(dataset)
        fprintf('data: %s\n', dataset{n}.name);
        
        %% read audio data with downsampling
        [x, fs, ~] = audioread_wrapper(dataset{n}.path, fs_d);
        
        %% cut the target segment
        x_t = cut(x, fs, dataset{n}.begin, dataset{n}.end);
        
       %% remove DC
        x_t = x_t - mean(x_t);
        
        %% STFT
        [V, ~, T_f] = spectrogram(x_t, w, M, N, fs);
        V = [V; conj(flipud(V(2:end-1, :)))];
        F_ft = (0:(N-1)).*(fs/N);
    
        figure;
        surf(T_f, F_ft(2:(N/2+1)), 2.*log(abs(V(2:(N/2+1), :))), 'edgecolor', 'none');
            axis tight;
            view(0, 90);
            set(gca, 'YScale', 'log');
            grid off;
            ylim([dataset{n}.ymin dataset{n}.ymax]);
        
        figure;
        surf(T_f, F_ft(2:(N/2+1)), 2.*log(abs(V(2:(N/2+1), :))), 'edgecolor', 'none');
            axis tight;
            view(0, 90);
            set(gca, 'YScale', 'log');
            grid off;
            
        %% CWT
        [s, F_wt] = linscaleset(N, fs, w0);
        omega = angfreqset(x_t, fs);
        g = morlet.scaledspectrum(omega, w0, s, 1/fs);
        W = morlet.transform(x_t, g);
        W = W';

        %% decimation
        dn = round(fs*dt);
        W = W(:, 1:dn:end);

        T_w = dataset{n}.begin + (0:dn:(length(x_t) - 1))./fs;
        
        figure;
        surf(T_w, F_wt, log(abs(W).^2), 'edgecolor', 'none');
            axis tight;
            view(0, 90);
            set(gca, 'YScale', 'log');
            grid off;
            ylim([dataset{n}.ymin dataset{n}.ymax]);
        
        figure;
        surf(T_w, F_wt, log(abs(W).^2), 'edgecolor', 'none');
            axis tight;
            view(0, 90);
            set(gca, 'YScale', 'log');
            grid off;
    end
end

%{
figure;
plot(F_ft(2:(N/2+1)) - F_wt);

t = randi(length(T_f));
[~, I] = min(abs(T_w - T_f(t)).^2);

maxval_ft = max(abs(V(2:(N/2+1), t)));
maxval_wt = max(abs(W(:, I)));

figure;
plot(F_ft(2:(N/2+1)), abs(V(2:(N/2+1), t))./maxval_ft); hold on;
plot(F_wt, abs(W(:, I))./maxval_wt, '-.m'); hold off;
xlim([dataset{n}.ymin dataset{n}.ymax]);
%}