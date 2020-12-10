function [x_d, fs, t] = audioread_wrapper(audiofilepath, fs_d)
    %% read audio file
    [x, fs] = audioread(audiofilepath);
    
    %% convert to mono if stereo
    if size(x, 2) == 2
        x_m = (x(:, 1) + x(:, 2))./2;
    else
        x_m = x;
    end
    
    %% down sampling
    x_d = x_m;
    
    while fs > fs_d
        x_d = x_d(1:2:end);
        fs = fs/2;
    end
    
    %% time axis
    t = (0:(length(x_d) - 1))./fs;
end