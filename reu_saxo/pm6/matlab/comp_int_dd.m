fs = 2760;

psd_path = '../results/bright1_03_8ms/cao_2/saxoplus_KL_phase_psd.fits';
psd_int_cao = fitsread(psd_path);

psd_path = '../results/bright1_03_8ms/cao_dd/saxoplus_KL_phase_psd.fits';
psd_dd_cao = fitsread(psd_path);

psd_path = '../results/bright1_03_8ms/dcao/saxoplus_KL_phase_psd.fits';
psd_int_dcao = fitsread(psd_path);

psd_path = '../results/bright1_03_8ms/dcao_dd/saxoplus_KL_phase_psd.fits';
psd_dd_dcao = fitsread(psd_path);

f = fitsread('../results/bright1_03_8ms/dcao/freq.fits');

%%

mode = 1;
figure()
% subplot(2,2,1)
semilogx(f,10*log10(psd_int_cao(:,mode)))
hold on;
semilogx(f,10*log10(psd_dd_cao(:,mode)))
semilogx(f,10*log10(psd_int_dcao(:,mode)))
semilogx(f,10*log10(psd_dd_dcao(:,mode)))
title('saxo+ standalone 1st KL residual PSD bright 1')
legend('int cao','dd cao','int dcao','dd dcao')
% legend_stand_g3 = sprintf('standalone: gain = 0.3 rms = %.2f',rms_standalone_g3);
% legend_dcao_g3 = sprintf('dcao: gain = 0.3 rms = %.2f',rms_dcao_g3);
% legend(legend_stand_g3,legend_dcao_g3)
% title('Tip residual, bright case, seing = 0.3", t0 = 2ms')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

%%
res_path = '../results/bright1_03_8ms/dcao/saxoplus_KL_res.fits';
res_int = fitsread(res_path);

res_path = '../results/bright1_03_8ms/dcao_dd/saxoplus_KL_res.fits';
res_dd = fitsread(res_path);
%%
figure()

plot(res_int(:,mode))
hold on;
plot(res_dd(:,mode))

title('saxo+ standalone 1st KL residual bright 1')
legend('int','dd')
rms(res_int(:,mode))
rms(res_dd(:,mode))