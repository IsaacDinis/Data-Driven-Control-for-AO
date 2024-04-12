fs = 1200;
TT_vibr = load('data/TT_pol_SPHERE.mat');
tilt_s1 = TT_vibr.TT_psol_s1(:,1);

tilt_s1_filt = bandpass(tilt_s1,[10,160],fs);
size_fft = 200;

[tilt_s1_psd, f] = compute_psd_welch(tilt_s1,size_fft,fs);
[tilt_s1_filt_psd, f] = compute_psd_welch(tilt_s1_filt,size_fft,fs);

t = 0:1/fs:length(tilt_s1)/fs-1/fs;

figure()
subplot(1,2,1)
plot(t,tilt_s1)
hold on
plot(t,tilt_s1_filt)
legend('raw tilt','filtered tilt')
xlabel('time (s)')
ylabel('tilt (mas)')
title('time signal')
make_it_nicer()
subplot(1,2,2)
semilogx(f,10*log10(tilt_s1_psd))
hold on;
semilogx(f,10*log10(tilt_s1_filt_psd))
legend('raw tilt','filtered tilt')
xlabel('freq (Hz)')
ylabel('tilt mag (dB)')
title('PSD')
make_it_nicer()
sgtitle('Bandpass tilt signal to only keep vibration')
%%
fs_sim = 3000;
tilt_s1_filt_sim = resample(tilt_s1_filt,fs_sim,fs);
mas2rad = 4.84814e-9;
rad2um = 1.65/2/pi;

%%
tilt_ol = fitsread('tilt_ol.fits');
tilt_s1_filt_sim = tilt_s1_filt_sim(1:length(tilt_ol));
phase2nm = 2.873136150686732;
factor = 4.848E-9*8E6;
full_tilt = tilt_ol+tilt_s1_filt_sim*factor*4;
[full_tilt_psd, f] = compute_psd_welch(full_tilt,size_fft,fs_sim);
[tilt_ol_psd, f] = compute_psd_welch(tilt_ol,size_fft,fs_sim);
figure()
semilogx(f,10*log10(full_tilt_psd))
hold on;
semilogx(f,10*log10(tilt_ol_psd))

fitswrite(tilt_s1_filt_sim*factor,'tilt_vibration.fits')

% to rad (tanj) *8m 