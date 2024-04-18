fs = 1200;
% TT_vibr = load('data/TT_pol_SPHERE.mat');
% tilt_s1 = TT_vibr.TT_psol_s1(:,1);
% tilt_s1 = wdenoise(double(tilt_s1));
% tilt_s1 = detrend(tilt_s1);
% tilt_s1 = tilt_s1(1:10000);
% tilt_s1 = tilt_s1-mean(tilt_s1);
% mas2um = 4.84e-9*8e6;
% tilt_s1 = tilt_s1*mas2um;
tilt_s1 = fitsread('tilt_vibration.fits');
tilt_s1 = tilt_s1(1:20000);
tilt_s1 = wdenoise(double(tilt_s1));
u = wgn(length(tilt_s1),1,0.5);
m = 8000;
N = length(tilt_s1);

T = [0:1:N-1]/fs;
toe = toeplitz(u(1:(end-1)));
PHI = tril(toe);
PHI = PHI(:,1:m);

theta = (PHI'*PHI)\ PHI' * tilt_s1(2:end);
y_hat = PHI*theta;
loss = sum((tilt_s1(2:end)-y_hat).^2);
fprintf(" Loss function of our prediction : %f \n",loss);
figure('Name','Estimated output','NumberTitle','off');
plot(T,tilt_s1,'DisplayName', 'measured');
hold on;
plot(T(2:end),y_hat, 'DisplayName', 'esitimated');
title('Estimated vs measured output');
xlabel('Time [s]');
ylabel('Amplitutde')
legend


size_fft = 500;

[tilt_s1_psd, f] = compute_psd_welch(tilt_s1,size_fft,fs);
[yhat_psd, f] = compute_psd_welch(y_hat,size_fft,fs);

figure()
semilogx(f,10*log10(tilt_s1_psd))
hold on;
semilogx(f,10*log10(yhat_psd))
legend('raw tilt','filtered tilt')
xlabel('freq (Hz)')
ylabel('tilt mag (dB)')

%%
% u = wgn(length(tilt_s1),1,0.5);
t = 0:1/fs:length(u)/fs-1/fs;
% u_test = wgn(length(tilt_s1),1,0.5);
u_test = wgn(50000,1,0.5);
yhat2 = conv(theta,u_test);
% yhat2 = yhat2(m:length(tilt_s1));
yhat2 = yhat2(m:length(yhat2)-m);

figure()
plot(tilt_s1)
hold on
plot(yhat2)
legend('raw tilt','filtered tilt')

[tilt_s1_psd, f] = compute_psd_welch(tilt_s1,size_fft,fs);
[yhat_psd, f] = compute_psd_welch(yhat2,size_fft,fs);


figure()
semilogx(f,10*log10(tilt_s1_psd))
hold on;
semilogx(f,10*log10(yhat_psd))
legend('initial tilt vibration','generated tilt vibration')
xlabel('freq (Hz)')
ylabel('PSD mag (dB)')
make_it_nicer()
fitswrite(yhat2,'tilt_vibration_gen.fits')

%%

% x = log10(f);
% y = 10*log10(tilt_s1_psd);
% 
% 
% 
% 
% 
% 
% 
% % curveFitter(x,y)
% [fitresult, gof] = create_fit_pol4(x, y');
% y_fit = feval(fitresult,x);
% 
% figure()
% plot(x,y)
% hold on
% plot(x,y_fit)