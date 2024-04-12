fs = 1200;
TT_vibr = load('data/TT_pol_SPHERE.mat');
tilt_s1 = TT_vibr.TT_psol_s1(:,1);
tilt_s1 = wdenoise(double(tilt_s1));
tilt_s1 = detrend(tilt_s1);
tilt_s1 = tilt_s1(1:5000);
tilt_s1 = tilt_s1-mean(tilt_s1)
mas2um = 4.84e-9*8e6;
tilt_s1 = tilt_s1*mas2um;
u = wgn(length(tilt_s1),1,0.5);
m = 1000;
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


size_fft = 200;

[tilt_s1_psd, f] = compute_psd_welch(tilt_s1,size_fft,fs);
[yhat_psd, f] = compute_psd_welch(y_hat,size_fft,fs);

figure()
semilogx(f,10*log10(tilt_s1_psd))
hold on;
semilogx(f,10*log10(yhat_psd))
legend('raw tilt','filtered tilt')
xlabel('freq (Hz)')
ylabel('tilt mag (dB)')

%% ARX
m = 20;
n = 20;
y = tilt_s1;
N = length(y);
u = wgn(length(tilt_s1),1,0.5);
toe_u = toeplitz(u(1:(end-1)));
PHI_u = tril(toe_u);
PHI_u = PHI_u(:,1:m);
toe_y = toeplitz(y(1:(end-1)));
PHI_y = tril(toe_y);
PHI_y = PHI_y(:,1:m);
PHI = [-PHI_y,PHI_u];

fprintf("Estimated parameters : \n");
theta = (PHI'*PHI)\ PHI' * y(2:end);
disp(theta);
y_hat = PHI*theta;
loss = sum((y(2:end)-y_hat).^2);
fprintf(" Loss function of our prediction : %f \n",loss);
figure('Name','Estimated output','NumberTitle','off');
plot(y(2:end),'DisplayName','measured');
hold on;
plot(y_hat,'DisplayName','estimated');
title('Estimated vs measured output');
xlabel('Time [ms]');
ylabel('Amplitutde')
legend
t = 0:1/fs:length(u)/fs-1/fs;
sys = tf(theta(21:40)',theta(1:20)',1/fs);
[yhat2,~] = lsim(sys,u,t);

figure()
plot(y)
hold on
plot(yhat2)

[y_fft, f] = compute_psd_welch(y,size_fft,fs);
[yhat2_fft, f] = compute_psd_welch(yhat2,size_fft,fs);
figure()
semilogx(f,20*log10(y_fft));
hold on;
semilogx(f,yhat2_fft)

%%
sys = ar(y,50,'ls','Ts',1/fs);
figure()
spectrum(sys)

p = aryule(y,10);
figure()
spectrum(sys)
%%
fs = 1200;
TT_vibr = load('data/TT_pol_SPHERE.mat');
tilt_s1 = TT_vibr.TT_psol_s1(:,1);
% tilt_s1 = wdenoise(double(tilt_s1));
tilt_s1 = detrend(tilt_s1);
tilt_s1 = tilt_s1-mean(tilt_s1);
mas2um = 4.84e-9*8e6;
tilt_s1 = tilt_s1*mas2um;
y = tilt_s1;
size_fft = 500;
y_low = lowpass(y,80,fs);
[y_fft, f] = compute_psd_welch(y,size_fft,fs);
y_fft = y_fft/max(y_fft);
f(end) = fs/2;
[b,a] = yulewalk(20,f'/fs*2,y_fft);
sys = tf(b,a,1/fs);

[mag,phase,wout] = bode(sys);

figure()
semilogx(f,20*log10(y_fft));
hold on;
semilogx(wout/2/pi,20*log10(squeeze(mag)))