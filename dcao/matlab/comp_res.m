tilt_saxo_standalone = fitsread('../data/tilt_saxo_standalone.fits')*1000;
tilt_saxoplus_standalone = fitsread('../data/tilt_saxoplus_standalone.fits')*1000;
res_saxo_mode_0 = fitsread('../data/res_saxo_mode_0.fits')*1000;
res_saxoplus_mode_0 = fitsread('../data/res_saxoplus_mode_0.fits')*1000;
% tilt_saxo_dcao = fitsread('../data/tilt_saxo_dcao.fits')*1000;
% tilt_saxoplus_dcao = fitsread('../data/tilt_saxoplus_dcao.fits')*1000;
tilt_saxo_dd = fitsread('../data/tilt_saxo_dd3.fits')*1000;
tilt_saxoplus_dd = fitsread('../data/tilt_saxoplus_dd3.fits')*1000;
T = 1;
fs = 1/0.00036231884057971015;
t = 0:1/fs:T-1/fs;
%%
figure()
plot(t,tilt_saxoplus_dd)
hold on;
plot(t,tilt_saxoplus_dcao)
%%
figure()
plot(t,tilt_saxoplus_dcao)
hold on;
plot(t,tilt_saxoplus_standalone)

%%
n = 20;
w = 100;

[psd_saxo_standalone,f] = compute_psd((tilt_saxo_standalone)',n,w,fs);
[psd_saxoplus_standalone,f] = compute_psd((tilt_saxoplus_standalone)',n,w,fs);
% [psd_saxo_dcao,f] = compute_psd((tilt_saxo_dcao)',n,w,fs);
% [psd_saxoplus_dcao,f] = compute_psd((tilt_saxoplus_dcao)',n,w,fs);
[psd_saxo_dd,f] = compute_psd((tilt_saxo_dd)',n,w,fs);
[psd_saxoplus_dd,f] = compute_psd((tilt_saxoplus_dd)',n,w,fs);

figure()
plot(f(1:end),psd_saxo_standalone(1:end))


figure()
plot(f(1:end),psd_saxoplus_dd(1:end))
hold on;
plot(f(1:end),psd_saxoplus_standalone(1:end))
title('tilt residual PSD, fs second stage = 2760')
legend('1st stage fs = 1380','1st stage fs = 690')
sum(psd_saxoplus_dd)
sum(psd_saxoplus_standalone)
%%
n = 20;
w = 100;

[psd_res_saxo_mode_0,f] = compute_psd((res_saxo_mode_0)',n,w,fs);
[psd_saxo_standalone,f] = compute_psd((tilt_saxo_standalone)',n,w,fs);

figure()
semilogx(f,10*log10(psd_res_saxo_mode_0))
hold on
semilogx(f,10*log10(psd_saxo_standalone))

%%
n = 20;
w = 100;

[psd_res_saxoplus_mode_0,f] = compute_psd((res_saxoplus_mode_0)',n,w,fs);
[psd_saxoplus_standalone,f] = compute_psd((tilt_saxoplus_standalone)',n,w,fs);

figure()
semilogx(f,10*log10(psd_res_saxoplus_mode_0))
hold on
semilogx(f,10*log10(psd_saxoplus_standalone))