function [psd, f] =  compute_psd(data,n_average,window_size,fs)
    n_modes = size(data,2);
    hann_window = hann(2*window_size+1);
    psd = zeros(2*window_size+1,n_modes);

    for mode = 1:n_modes
        for i = 1:n_average
            data_xcorr = xcorr(data(1+(i-1)*window_size:1+i*window_size,mode),'unbiased');
            data_xcorr_w = hann_window.*data_xcorr;
            psd(:,mode) = psd(:,mode) + abs(fft(data_xcorr_w));
        end
    end
    psd = psd/n_average;
    f = fs/(2*window_size+1)*(0:2*window_size);
    
    % remove frequencies above nyquist
    f = f(1:window_size+1);
    psd = psd(1:window_size+1,:);
end