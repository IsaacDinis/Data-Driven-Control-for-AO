function [psd_average, f, psd_std] =  compute_psd2(data,n_average,window_size,fs)
    n_modes = size(data,2);
    hann_window = hann(2*window_size+1);
    psd = zeros(2*window_size+1,n_modes,n_average);

    for mode = 1:n_modes
        for i = 1:n_average
            data_xcorr = xcorr(data(1+(i-1)*window_size:1+i*window_size,mode),'biased');
            data_xcorr_w = hann_window.*data_xcorr;
            psd(:,mode,i) = abs(fft(data_xcorr_w)/window_size);
        end
    end
    psd_average = sum(psd,3)/n_average;
    psd_std = std(psd,0,3);
    f = fs/(2*window_size+1)*(0:2*window_size);
    
    % remove frequencies above nyquist
    f = f(1:window_size+1);
    psd_average = psd_average(1:window_size+1,:);
    psd_std = psd_std(1:window_size+1,:);
end