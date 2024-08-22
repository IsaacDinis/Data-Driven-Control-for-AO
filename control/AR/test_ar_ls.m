fs = 4000;


tilt_res = fitsread('tilt_res.fits');
tilt_applied = fitsread('tilt_applied.fits');
x = tilt_applied(1:end-2)+tilt_res(1+2:end);
x = x(1:10000)';

m = 5;
N = length(x);

T = [0:1:N-1]/fs;
toe = toeplitz(x(1:end-1));
PHI = tril(toe);
PHI = PHI(:,1:m);

theta = (PHI'*PHI)\ PHI' * x(2:end);
y_hat = PHI*theta;
loss = sum((x(2:end)-y_hat).^2);
fprintf(" Loss function of our prediction : %f \n",loss);

figure('Name','Estimated output','NumberTitle','off');
plot(T(2:end),x(2:end),'DisplayName', 'true');
hold on;
plot(T(2:end),y_hat, 'DisplayName', 'prediction');
plot(T(2:end),x(1:end-1), 'DisplayName', 'measured');
title('Estimated vs measured output');
xlabel('Time [s]');
ylabel('Amplitutde')
legend

% 
% size_fft = 500;
% 
% [tilt_s1_psd, f] = compute_psd_welch(tilt_s1,size_fft,fs);
% [yhat_psd, f] = compute_psd_welch(y_hat,size_fft,fs);
% 
% figure()
% semilogx(f,10*log10(tilt_s1_psd))
% hold on;
% semilogx(f,10*log10(yhat_psd))
% legend('raw tilt','filtered tilt')
% xlabel('freq (Hz)')
% ylabel('tilt mag (dB)')
