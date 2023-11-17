function [psd, f] =  compute_psd_fft(data,n_average,window_size,fs)
    n_modes = size(data,2);
    psd = zeros(2*window_size+2,n_modes);

    for mode = 1:n_modes
        for i = 1:n_average
            data_w = data(1+(i-1)*window_size:1+i*window_size,mode);
            data_w = [data_w;zeros(size(data_w))];
            psd(:,mode) = psd(:,mode) + abs(fft(data_w)/size(data_w,1)).^2;
        end
    end
    psd = psd/n_average;
    f = fs/(2*window_size+1)*(0:2*window_size);
    
    % remove frequencies above nyquist
    f = f(1:window_size+1);
    psd = psd(1:window_size+1,:);
end