%{
N = 5;

gam_0 = 2 * rand(N, 1);
gam = 2 * rand(N, 2);

disp([gam_0 gam]);
main(gam_0, gam);
%}

%{
>> SPEのベースライン
gam_0 = 0.76;
gam = [1.14 1.00];

>> 最初にランダムで見つけたやつ
gam_0 = [0.76; 0.510];
gam = [1.14 1.00; 1.094 0.509];

>> 良くも悪くも普通. キレイになるが消しすぎる.
gam_0 = 1.0;
gam = [1.0 1.0 1.0 1.0];

>> ノイジーだが内声は残る
gam_0 = 1.0;
gam = [1.0 0.1 1.0 1.0];

>> 上記がちょっと改善
gam_0 = 1.0;
gam = [1.0 0.1 1.0 1.9];

>> 上記が更に改善するが内声が消え始める
gam_0 = 1.0;
gam = [1.0 0.1 1.9 1.9];

>> 辛うじて内声が保持できる水準
gam_0 = 0.45;
gam = [0.45 0.45 0.45 0.45];

>> ほぼ消える
gam_0 = 1.5;
gam = [1.5 1.5 1.5 1.5];

>> もう少しノイジーにする？
gam_0 = 0.3;
gam = [0.3 1.0 1.0 0.5 0.5];

>> そしてこれが最終的な調整結果
gam_0 = 0.1;
gam = [0.3 1.0 1.0 0.5 0.5];

>> これだとちょっとだけノイジーすぎるかも
gam_0 = 0.05;
gam = [0.3 1.0 1.0 0.5 0.5];

>> 合体
gam_0 = [0.6; 0.7];
gam = [0.4 1.9 1.0 0.1; 0.3 1.8 1.0 0.1];

>> 合体 2
gam_0 = [1.0; 0.45];
gam = [1.0 0.1 1.0 1.0; 0.45 0.1 0.45 0.45];
%}

function main(gam_0, gam)
    %% parameters
    if nargin ~= 2
        gam_0 = [1.0; 0.45];
        gam = [1.0 0.1 1.0 1.0; 0.45 0.1 0.45 0.45];
    end
    
    fs_d = 25000;
    
    N = 4096;
    
    %{
    M = round(0.98*N);
    w = gausswindow(N);
    %}
    
    %%{
    w0 = 64;
    dt = 0.005;
    %}
    
    g = @(x, gam) x.^gam;
    
    reffreq = 440;
    numsemi = 12;
    
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
        
       %% reference pitch search
        reffreq = refpitch(x_t, fs, reffreq, numsemi);
        fprintf('A4 = %3.3f Hz\n', reffreq);
        
       %% STFT
        %{
        [V, ~, T_f] = spectrogram(x_t, w, M, N, fs);
        V = [V; conj(flipud(V(2:end-1, :)))];
        %}

       %% CWT
        %%{
        [s, ~] = linscaleset(N, fs, w0);
        omega = angfreqset(x_t, fs);
        gm = morlet.scaledspectrum(omega, w0, s, 1/fs);
        V = morlet.transform(x_t, gm);

        dn = round(fs*dt);
        V = V(1:dn:end, :);

        V = V';

        V = [1e-10.*ones(1, size(V, 2)); V];
        V = [V; conj(flipud(V(2:end-1, :)))];

        T_f = (0:dn:(length(x_t) - 1))./fs;
        %}

        %% dimension
        F = (0:(N-1)).*(fs/N);
        
        %{
        %% Wavelet-based MLC
        T = 1./F(2:(N/2 + 1));
        s = T.*(w0 + sqrt(2 + w0^2))./(4*pi);
        omega = angfreqset(ones(N, 1), fs);
        gm = morlet.scaledspectrum(omega, w0, s, 1/fs);
        
        t = 160;
        C = morlet.transform(V(:, t), gm);
        %C(C < 0) = 0;
        C(:, F(2:(N/2 + 1)) > q_k) = 0;
        W = log(abs(V(2:(N/2 + 1), t))) .* log(abs(C(1, :)'));
        
        figure;
        subplot(5, 1, 1); plot(F(2:(N/2 + 1)), abs(V(2:(N/2 + 1), t)));
        title(num2str(T_f(t) + dataset{n}.begin));
        subplot(5, 1, 2); plot(abs(C(1, :)'));
        subplot(5, 1, 3); plot(F(2:(N/2 + 1)), abs(W));
        subplot(5, 1, 4); plot(F(2:(N/2 + 1)), abs(V(2:(N/2 + 1), t))); xlim([0 1000]);
        subplot(5, 1, 5); plot(F(2:(N/2 + 1)), abs(W)); xlim([0 1000]);
        
        figure;
        surf(T_f(t-50:t+50) + dataset{n}.begin, F(1:(N/2+1)), abs(V(1:(N/2+1), t-50:t+50)), 'edgecolor', 'none');
            title(sprintf('%s (ensemble, fs = %1.2f)', dataset{n}.name, fs));
            axis tight;
            view(0, 90);
            grid off;
            ylim([0 1000]);
        
        figure;
        surf(1:4096, 1:2048, abs(C'), 'edgecolor', 'none');
        axis tight;
        view(0, 90);
        %}
        
        %% energy freq lower-bound
        [f_k, q_k, S] = frangepower(x_t, fs, w0, dt);
        fprintf('99%% frequency energy cutoff line: %3.3f Hz / %3.3f Hz\n', f_k, q_k);
        
        %% MLC
        W = cell(1, length(gam_0));
        W_q = cell(1, length(gam_0));
        
        for m=1:length(gam_0)
            [W_m, z, z_q] = multilayeredceps(F, f_k, q_k, V, gam_0(m, :), gam(m, :), g);
            
            [z, F_q] = freqquant(z, fs, f_k, q_k, reffreq, numsemi);
            z_q = quefquant(z_q, fs, f_k, q_k, reffreq, numsemi);
            
            W{m} = W_m;
            W_q{m} = z .* z_q;
        end
        
        W_e = 0;
        E = sum(abs(V(:).^2))./size(V, 2);
        for m=1:length(gam_0)
            E_m = sum(abs(W{m}(:)).^2)./size(W{m}, 2);
            
            W{m} = W{m} .* sqrt(E/E_m);
            
            W_e = W_e + 2.*log(abs(W{m}));
        end
        
        %% plot
        numfig = 2;
        subfignum = numfig*(n - 1);

        T_f = T_f + dataset{n}.begin;
        plabel = pitchlabel(F_q, reffreq);
        
        figure;
        surf(T_f, F(1:(N/2+1)), 2.*log(abs(V(1:(N/2+1), :))), 'edgecolor', 'none');
            title(sprintf('%s (ensemble, fs = %1.2f)', dataset{n}.name, fs));
            axis tight;
            view(0, 90);
            set(gca, 'YScale', 'log');
            grid off;
            
        
        for m=1:length(gam_0)
            W_m = W{m};
            W_m = 2.*log(abs(W_m));
            W_m(W_m == -Inf) = -100;
            
            figure;
            surf(T_f, F(1:(N/2+1)), W_m(1:(N/2+1), :), 'edgecolor', 'none');
                title(sprintf('%s (m = %d, fs = %1.2f)', dataset{n}.name, m, fs));
                axis tight;
                view(0, 90);
                set(gca, 'YScale', 'log');
                yticks(F_q);
                grid off;
                hold on
            
            surf([T_f(1) T_f(end)], F_q, ones(length(F_q), 2),...
                'edgecolor', [0.95 0.8 0.95], 'linestyle', '--', 'facealpha', 0);
                grid off
                view(0, 90);
                axis tight;
                set(gca, 'YScale', 'log');
                hold off
                
                yticklabels(plabel);
                ylim([f_k dataset{n}.ymax]);
            
            figure;
            surf(T_f, F(1:(N/2+1)), W_m(1:(N/2+1), :), 'edgecolor', 'none');
                title(sprintf('%s (m = %d, fs = %1.2f)', dataset{n}.name, m, fs));
                axis tight;
                view(0, 90);
                set(gca, 'YScale', 'log');
                grid off;
        end
       
       %% plot
        W_t = W_e;
        W_t(W_t == -Inf) = -100;
        
        figure;
        surf(T_f, F(1:(N/2+1)), W_t(1:(N/2+1), :), 'edgecolor', 'none');
            title(sprintf('%s (ensemble, fs = %1.2f)', dataset{n}.name, fs));
            axis tight;
            view(0, 90);
            set(gca, 'YScale', 'log');
            grid off;
        
        figure;
        surf(T_f, F(1:(N/2+1)), W_t(1:(N/2+1), :), 'edgecolor', 'none');
            title(sprintf('%s (ensemble, fs = %1.2f)', dataset{n}.name, fs));
            axis tight;
            view(0, 90);
            set(gca, 'YScale', 'log');
            yticks(F_q);
            grid off;
            hold on;
           
        surf([T_f(1) T_f(end)], F_q, ones(length(F_q), 2),...
            'edgecolor', [0.95 0.8 0.95], 'linestyle', '--', 'facealpha', 0);
            grid off
            view(0, 90);
            axis tight;
            set(gca, 'YScale', 'log');
            hold off
                
            yticklabels(plabel);
            ylim([f_k dataset{n}.ymax]);
        
       %% post-processing
        W_p = exp(W_e);
        
        % silence part masking
        W_p = bsxfun(@times, W_p, S');
        
        % peak picking verifying with original scalogram
        eps_c = 40;
        
        W_t = W_p .* 0;
        V_p = abs(V).^2;
        
        I_cutoff = max(find((F - q_k) < 0));
        L_cutoff = min(find((F - f_k) > 0));
        
        for t=1:size(W_t, 2)
            %E_w = log(W_p(:, t));
            %E_w(E_w == -Inf) = -20;
            E_w = W_p(:, t);
            
            %E_v = log(V_p(:, t));
            %E_v(E_v == -Inf) = -20;
            E_v = V_p(:, t);
            
            [~, locs_w] = findpeaks(E_w(L_cutoff:I_cutoff));
            [~, locs_v] = findpeaks(E_v(L_cutoff:I_cutoff));
            
            C_w = 1200*log2(F(locs_w + L_cutoff - 1)./440);
            C_v = 1200*log2(F(locs_v + L_cutoff - 1)./440);
            
            I = ones(length(C_w), 1) == 0;
            for i=1:length(I)
                if min(abs(C_w(i) - C_v)) < eps_c
                    I(i) = 1;
                end
            end
            
            W_t(locs_w(I) + L_cutoff - 1, t) = W_p(locs_w(I) + L_cutoff - 1, t);
        end
        
        %{
        figure;
        subplot(2, 1, 1);
        plot(F, E_w); hold on;
        stem(F(locs_w(I) + L_cutoff - 1), E_w(locs_w(I) + L_cutoff - 1), 'color', 'm'); hold off;
        xlim([F(L_cutoff) F(I_cutoff)]);
        
        subplot(2, 1, 2);
        plot(F, E_v); hold on;
        
        E_c = dct(E_v);
        E_c(round(N/5):end) = 0;
        E_c = idct(E_c);
        
        plot(F, E_c, '-.m'); hold off;
        
        xlim([F(L_cutoff) F(I_cutoff)]);
        title(num2str(T_f(t)));
        %}
        
        %{
        for t=1:size(W_p, 2)
            [~, locs] = findpeaks(W_p(L_cutoff:I_cutoff, t));
            locs = locs + L_cutoff - 1;
            
            W_t(locs, t) = W_p(locs, t);
        end
        %}
        
        W_p = W_t;
        
        % low dB cutoff
        W_p(log(W_p) <= -20) = 0;
        
        % eliminate quasi-F0
        W_p = W_p .* abs(V).^2;
        
        % quantisation
        [W_q, F_wq] = freqquant(W_p, fs, f_k, ceil(q_k * 2^(300/1200)), 440, 120);
        
        % time-smoothing filtering
        %{
        H = ones(1, 3);
        H = H./length(H);
        W_q = filter2(H, W_q);
        %}
        
        % vertical thresholding (masking clustering)
        eps_l = log(1/400);
        eps_i = 9;
        [~, I] = max(W_q, [], 1);
        W_t = W_q .* 0;
        for t=1:size(W_t, 2)
            E = W_q(:, t);
            
            r = E./E(I(t));
            E(log(r) < eps_l) = 0;
            
            I_c = find(E ~= 0);
            
            if ~isempty(I_c)
                c = I_c .* 0;
                K = 1;
                c(1) = K;

                for i=2:length(c)
                    if (I_c(i) - I_c(i - 1)) < eps_i
                        c(i) = c(i - 1);
                    else
                        K = K + 1;
                        c(i) = K;
                    end
                end

                for k=1:K
                    I_k = I_c(c == k);

                    [P, I_p] = max(E(I_k));
                    I_p = I_k(I_p);

                    E(I_k) = 0;
                    E(I_p) = P;
                end
            end
            
            W_t(:, t) = E;
        end
        W_q = W_t;
        
         %% plot
        W_t = log(W_p);
        W_t(W_t == -Inf) = -100;
        
        figure;
        surf(T_f, F(1:(N/2+1)), W_t(1:(N/2+1), :), 'edgecolor', 'none');
            title(sprintf('%s (ensemble, fs = %1.2f)', dataset{n}.name, fs));
            axis tight;
            view(0, 90);
            set(gca, 'YScale', 'log');
            grid off;
        
        figure;
        surf(T_f, F(1:(N/2+1)), W_t(1:(N/2+1), :), 'edgecolor', 'none');
            title(sprintf('%s (ensemble, fs = %1.2f)', dataset{n}.name, fs));
            axis tight;
            view(0, 90);
            set(gca, 'YScale', 'log');
            yticks(F_q);
            grid off;
            hold on;
           
        surf([T_f(1) T_f(end)], F_q, ones(length(F_q), 2),...
            'edgecolor', [0.3 0.3 0.85], 'linestyle', '--', 'facealpha', 0);
            grid off
            view(0, 90);
            axis tight;
            set(gca, 'YScale', 'log');
            hold off
                
            yticklabels(plabel);
            ylim([f_k dataset{n}.ymax]);
        
       %% plot
        W_t = log(W_q);
        W_t(W_t == -Inf) = -100;
        
        figure;
        surf(T_f, F_wq, W_t, 'edgecolor', 'none');
            title(sprintf('%s (ensemble, fs = %1.2f)', dataset{n}.name, fs));
            axis tight;
            view(0, 90);
            set(gca, 'YScale', 'log');
            yticks(F_q);
            grid off;
            hold on;
           
        surf([T_f(1) T_f(end)], F_q, ones(length(F_q), 2),...
            'edgecolor', [0.3 0.3 0.85], 'linestyle', '--', 'facealpha', 0);
            grid off
            view(0, 90);
            axis tight;
            set(gca, 'YScale', 'log');
            hold off
                
            yticklabels(plabel);
            ylim([f_k dataset{n}.ymax]);
            
        %% pitch tracking
        [row, col] = ind2sub(size(W_q), find(W_q(:) ~= 0));
        
        p = zeros(2, length(row));
        col_I = unique(col);
        for i=1:length(col_I)
            p(1, col == col_I(i)) = T_f(col_I(i));
        end
        
        I = unique(row);
        for i=1:length(I)
            c = 1200 * log2(F_wq(I(i)) / 440);
            p(2, row == I(i)) = c;
        end
        
        [~, I] = sort(p(1, :));
        p(1, :) = p(1, I);
        p(2, :) = p(2, I);
        
        eps_t = 0.1;
        eps_c = 100;
        c = pitchtracking(p, W_q, eps_t, eps_c, col_I);
        
        eps_n = 80;
        
        K = length(unique(c));
        for k=1:K
            I = find(c == k);
            
            if length(I) <= eps_n
                c(I) = 0;
            end
        end
        
        I = c == 0;
        K = length(unique(c)) - 1;
        
        figure;
        gscatter(p(1, ~I), p(2, ~I), c(~I)); title(sprintf('%d clusters\n', K));
        legend('off');
        
        %{
        I = c ~= 0;
        kset = unique(c(I));
        figure; 
        for k=1:K
            I = find(c == kset(k));
            
            subplot(2, 1, 1);
            gscatter(p(1, :), p(2, :), c); hold on
            legend('off');
            scatter(p(1, I), p(2, I)); hold off
            title(sprintf('cluster %d (%d data points)', kset(k), length(I)));
            
            subplot(2, 1, 2); plot(p(1, I), p(2, I));
            
            pause;
        end
        %}
        
        %% merge overlap clustering (interpolation and distance)
        % extract time-domain range
        
        
        [row, col] = ind2sub(size(W_q), find(W_q(:) ~= 0));
        I = c ~= 0;
        K = length(unique(c(I)));
        E_c = zeros(2, K);
        recbound = zeros(4, K); %left, right, top, bottom
        kset = unique(c(I));
        W_t = W_q .* 0;
        
        for k=1:K
            c_k = kset(k);
            
            I_k = find(c == c_k);
            
            xy = zeros(2, length(I_k));
            
            for i=1:length(I_k)
                E_c(1, k) = c_k;
                E_c(2, k) = E_c(2, k) + W_q(row(I_k(i)), col(I_k(i)));
                
                xy(1, i) = col(I_k(i));
                xy(2, i) = row(I_k(i));
                
                W_t(row(I_k(i)), col(I_k(i))) = W_q(row(I_k(i)), col(I_k(i)));
            end
            
            recbound(1, k) = min(xy(1, :));
            recbound(2, k) = max(xy(1, :));
            recbound(3, k) = max(xy(2, :));
            recbound(4, k) = min(xy(2, :));
        end
        
        [~, I_E] = sort(E_c(2, :), 'desc');
        E_c(1, :) = E_c(1, I_E);
        E_c(2, :) = E_c(2, I_E);
        recbound(1, :) = recbound(1, I_E);
        recbound(2, :) = recbound(2, I_E);
        recbound(3, :) = recbound(3, I_E);
        recbound(4, :) = recbound(4, I_E);
        
        %{
        W_t = log(W_t);
        W_t(W_t == -Inf) = -100;
        
        figure;
        surf(T_f, F_wq, W_t, 'edgecolor', 'none');
            title(sprintf('%s (ensemble, fs = %1.2f)', dataset{n}.name, fs));
            axis tight;
            view(0, 90);
            set(gca, 'YScale', 'log');
            yticks(F_q);
            grid off;
        %}
    end
end