% tilt_saxo_standalone = fitsread('../data/tilt_saxo_standalone.fits')*1000;
tilt_saxoplus_standalone = fitsread('../data/tilt_saxoplus_standalone.fits')*1000;
tilt_saxoplus_dcao = fitsread('../data/tilt_saxoplus_dcao.fits')*1000;
% command_rms_dcao = fitsread('../data/command_0_dcao.fits')*1000;
% command_rms_standalone = fitsread('../data/command_0_standalone.fits')*1000;

res_0_saxoplus_standalone = fitsread('../data/res_0_saxoplus_standalone.fits')*1000;
res_0_saxoplus_dcao = fitsread('../data/res_0_saxoplus_dcao.fits')*1000;
command_0_dcao = fitsread('../data/command_0_dcao.fits')*1000;
command_0_standalone = fitsread('../data/command_0_standalone.fits')*1000;

T = 1;
fs = 1/0.00036231884057971015;
t = 0:1/fs:T-1/fs;
%%
% figure()
% plot(t,tilt_saxoplus_dd)
% hold on;
% plot(t,tilt_saxoplus_dcao)
% %%

%%

figure()
plot(command_0_dcao)
hold on;
plot(command_0_standalone)
%%
n = 20;
w = 100;

[psd_saxoplus_standalone,f] = compute_psd((res_0_saxoplus_standalone)',n,w,fs);
[psd_saxoplus_dcao,f] = compute_psd((res_0_saxoplus_dcao)',n,w,fs);



figure()
semilogx(f,10*log10(psd_saxoplus_dcao))
hold on
semilogx(f,10*log10(psd_saxoplus_standalone))

title('2nd stage phase KL 0 residual PSD')
legend('dcao','standalone','Interpreter','latex','location','northeast');
ylabel('Magnitude (dB)')
xlabel('Frequency (Hz)')
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
sum(psd_saxoplus_standalone)
sum(psd_saxoplus_dcao)

%%
n = 20;
w = 100;

[psd_saxoplus_standalone,f] = compute_psd((tilt_saxoplus_standalone)',n,w,fs);
[psd_saxoplus_dcao,f] = compute_psd((tilt_saxoplus_dcao)',n,w,fs);



figure()
semilogx(f,10*log10(psd_saxoplus_dcao))
hold on
semilogx(f,10*log10(psd_saxoplus_standalone))

title('2nd stage phase tilt residual PSD')
legend('dcao','standalone','Interpreter','latex','location','northeast');
ylabel('Magnitude (dB)')
xlabel('Frequency (Hz)')
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
sum(psd_saxoplus_standalone)
sum(psd_saxoplus_dcao)
%%
n = 20;
w = 100;

[psd_saxoplus_standalone,f] = compute_psd((command_0_standalone)',n,w,fs);
[psd_saxoplus_dcao,f] = compute_psd((command_0_dcao)',n,w,fs);



figure()
semilogx(f,10*log10(psd_saxoplus_dcao))
hold on
semilogx(f,10*log10(psd_saxoplus_standalone))

title('2nd stage command KL 0 PSD')
legend('dcao','standalone','Interpreter','latex','location','northeast');
ylabel('Magnitude (dB)')
xlabel('Frequency (Hz)')
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
% sum(psd_saxoplus_dd)
% sum(psd_saxoplus_standalone)
