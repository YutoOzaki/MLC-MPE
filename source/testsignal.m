function testsignal
    fs = 44100;
    
    a = 0.3;
    
    f0 = 215 + 10*rand;
    
    t = (0:(fs-1))./fs;
    
    y = [];
    w = hann(fs)';
    for i=0:12
        f = f0*2^(i/12);
        
        x = w.*a.*sin(2.*pi.*f.*t);
        
        y = [y x];
    end
    
    audiowrite('testsignal.wav', y, fs);
    fprintf('A3 = %3.3f\n', f0);
    
    figure; spectrogram(y, hann(8192), 8192*0.75, 8192, fs, 'yaxis'); ylim([(f0-10) (2*f0+10)]./1000);
end