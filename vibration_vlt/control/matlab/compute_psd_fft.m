function [psd, f, spectrogram] =  compute_psd_fft(data,fft_size,fs)
    if mod(fft_size,2) == 0
        fft_size = fft_size+1;
    end
    window_size = fft_size*2-1;
    n_modes = size(data,2);
    n_frames = floor(size(data,1)/window_size);
    spectrogram = zeros(fft_size,n_frames,n_modes);

    for mode = 1:n_modes
        for i = 1:n_frames
            data_w = data(1+(i-1)*window_size:i*window_size,mode);
            psd_w = abs(fft(data_w))/window_size;
            psd_w = psd_w(1:fft_size);
            psd_w(2:end) = 2*psd_w(2:end);
            spectrogram(:, i, mode) = psd_w;
        end
    end
    psd = squeeze(mean(spectrogram,2));
    f = fs*(0:(window_size/2))/window_size;
    
    % remove frequencies above nyquist
end

% Y = fft(X);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;