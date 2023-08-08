% tilt_saxo_standalone = fitsread('../data/tilt_saxo_standalone.fits')*1000;
tilt_saxoplus_standalone = fitsread('../data/tilt_saxoplus_dcao_offset.fits')*1000;
% res_saxo_mode_0 = fitsread('../data/res_saxo_mode_0.fits')*1000;
% res_saxoplus_mode_0 = fitsread('../data/res_saxoplus_mode_0.fits')*1000;
% tilt_saxo_dcao = fitsread('../data/tilt_saxo_dcao.fits')*1000;
% tilt_saxoplus_dcao = fitsread('../data/tilt_saxoplus_dcao.fits')*1000;
% tilt_saxo_dd = fitsread('../data/tilt_saxo_dd5.fits')*1000;
tilt_saxoplus_dcao = fitsread('../data/tilt_saxoplus_dcao4.fits')*1000;
T = 1;
fs = 1/0.00036231884057971015;
t = 0:1/fs:T-1/fs;
%%
% figure()
% plot(t,tilt_saxoplus_dd)
% hold on;
% plot(t,tilt_saxoplus_dcao)
% %%
% figure()
% plot(t,tilt_saxoplus_dcao)
% hold on;
% plot(t,tilt_saxoplus_standalone)

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
% sum(psd_saxoplus_dd)
% sum(psd_saxoplus_standalone)
%%
% n = 20;
% w = 100;
% 
% [psd_res_saxo_mode_0,f] = compute_psd((res_saxo_mode_0)',n,w,fs);
% [psd_saxo_standalone,f] = compute_psd((tilt_saxo_standalone)',n,w,fs);
% 
% psd_saxo_standalone = psd_saxo_standalone/max(psd_saxo_standalone);
% psd_res_saxo_mode_0 = psd_res_saxo_mode_0/max(psd_res_saxo_mode_0);
% 
% figure()
% semilogx(f,10*log10(psd_res_saxo_mode_0))
% hold on
% semilogx(f,10*log10(psd_saxo_standalone))
% 
% title('2nd stage tilt disturbance PSD')
% legend('reconstructed 1st KL from WFS','phase tilt','Interpreter','latex','location','northeast');
% ylabel('Magnitude (dB)')
% xlabel('Frequency (Hz)')
% make_it_nicer()
% set(gcf, 'Position',  [100, 100, 700, 450])
% set(gcf,'PaperType','A4')

%%
% n = 20;
% w = 100;
% 
% [psd_res_saxoplus_mode_0,f] = compute_psd((res_saxoplus_mode_0)',n,w,fs);
% [psd_saxoplus_standalone,f] = compute_psd((tilt_saxoplus_standalone)',n,w,fs);
% 
% psd_res_saxoplus_mode_0 = psd_res_saxoplus_mode_0/max(psd_res_saxoplus_mode_0);
% psd_saxoplus_standalone = psd_saxoplus_standalone/max(psd_saxoplus_standalone);
% 
% figure()
% semilogx(f,10*log10(psd_res_saxoplus_mode_0))
% hold on
% semilogx(f,10*log10(psd_saxoplus_standalone))
% 
% title('2nd stage tilt residual PSD')
% legend('reconstructed 1st KL from WFS','phase tilt','Interpreter','latex','location','southwest');
% ylabel('Magnitude (dB)')
% xlabel('Frequency (Hz)')
% make_it_nicer()
% set(gcf, 'Position',  [100, 100, 700, 450])
% set(gcf,'PaperType','A4')