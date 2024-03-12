function [psd, f, spectrogram] =  compute_psd_welch(data,fft_size,fs)
    n_modes = size(data,2);
    % if mod(window_size,2) ~= 0
    %     window_size = window_size+1;
    % end
    if mod(fft_size,2) == 0
        fft_size = fft_size+1;
    end
    window_size = fft_size*2-1;

    n_frames = floor(size(data,1)/window_size)*2-1;
    spectrogram = zeros(fft_size,n_frames,n_modes);

    window = hamming(window_size)*2;

    for mode = 1:n_modes
        for i = 1:n_frames
            data_w = data(1+(i-1)*fft_size:(i-1)*fft_size+window_size, mode);
            data_w = data_w.*window;
            psd_w = abs(fft(data_w))/window_size;
            psd_w = psd_w(1:fft_size);
            psd_w(2:end) = 2*psd_w(2:end);
            spectrogram(:, i, mode) = psd_w;
        end
    end

    psd = squeeze(mean(spectrogram,2));
    
    f = fs*(0:(window_size/2))/window_size;
    
end