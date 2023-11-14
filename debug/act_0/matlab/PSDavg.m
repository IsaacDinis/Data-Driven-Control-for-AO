function [PSD, freq, cumulativePSD] = PSDavg(time, datat, sectionLength)
    % Function [PSD, freq, cumulativePSD] = PSDavg(time, data, sectionLength)
    %--------------------------------------------------------------
    % Compute the PSD of the data set.
    % !!! Returns the squared modulus of the FFT of the signal. !!!
    %
    % INPUT
    %  time: acquisition time of the data set 
    %    /!\ time must be uniformly sampled /!\
    %  data: data set
    %  sectionLength: length of sections used for computing PSDs,
    %  before they are averaged
    %
    % OUTPUT
    %  - Signal PSD
    %  - Frequency bins of PSD
    %  - Cumulative PSD - cPSD(end) = var(datat)
    %
    % FIX ME
    %
    % -------------------------------------------------------------
    
    Ndata     = length(datat);
    Navg      = floor(Ndata/sectionLength);  % round to next smaller integer (remove last bunch of data)
    Navg      = max(Navg, 1);  % in case the section length does not allow for Navg = 2
    if Navg < 1.5
       sectionLength = Ndata; 
    end
    
    time      = time - time(1);
    N         = sectionLength;
    dt        = mean(diff(time));
    max_ntime = N*dt;
    freq      = linspace(1./max_ntime, 1/dt, sectionLength);
    
    PSDi = zeros(Navg,N);
    for i=1:Navg
        lim_inf = (i-1)*N+1;
        lim_sup = i*N;​
        ddd = datat(lim_inf:lim_sup);
        ddd = ddd - mean(ddd); % Remove mean of each section
        xdft      = fft(ddd);
        PSDi(i,:) = (1/N) * abs(xdft).^2;
    end
    
    PSD = mean(PSDi, 1);
    cumulativePSD = 2*cumsum(PSD(1,1:floor(N/2)))/N;
    
    PSD = PSD(1:N/2);
    freq = freq(1:N/2);
end % PSDavg


















