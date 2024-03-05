Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L)*T;        % Time vector

S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);%+0.3*sin(2*pi*100*t);
X = S  ;%+ 2*randn(size(t));

% plot(1000*t(1:50),X(1:50))
title("Signal Corrupted with Zero-Mean Random Noise")
xlabel("t (milliseconds)")
ylabel("X(t)")

Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure()
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")

% sqrt(sum(P1)/f(2))
% rms(X)*2
% rms(S)*sqrt(2)

[psd_fft, f_fft] =  compute_psd_fft(S',751,Fs);
hold on
plot(f_fft,psd_fft)


[psd_welch, f_welch] =  compute_psd_welch(S',751,Fs);
plot(f_welch,psd_welch)
legend('matlab','fft','welch')

sum(psd_welch)
sum(P1)
sum(psd_fft)