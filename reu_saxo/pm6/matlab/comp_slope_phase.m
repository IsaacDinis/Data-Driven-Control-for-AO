fs = 2760;

psd_path = '../results/bright1_03_8ms/cao/saxoplus_KL_psd.fits';
psd_slope = fitsread(psd_path);

psd_path = '../results/bright1_03_8ms/cao/saxoplus_KL_phase_psd.fits';
psd_phase = fitsread(psd_path);

f = fitsread('../results/bright1_03_8ms/dcao/freq.fits');

%%

mode = 2;
figure()
% subplot(2,2,1)
semilogx(10*log10(psd_slope(:,mode)))
hold on;
semilogx(10*log10(psd_phase(:,mode)))

title('saxo+ standalone 1st KL residual PSD bright 1')
legend('slopes','phase')
% legend_stand_g3 = sprintf('standalone: gain = 0.3 rms = %.2f',rms_standalone_g3);
% legend_dcao_g3 = sprintf('dcao: gain = 0.3 rms = %.2f',rms_dcao_g3);
% legend(legend_stand_g3,legend_dcao_g3)
% title('Tip residual, bright case, seing = 0.3", t0 = 2ms')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

%%
res_path = '../results/bright1_03_8ms/cao/saxoplus_KL_res.fits';
res_slope = fitsread(res_path);

res_path = '../results/bright1_03_8ms/cao/saxoplus_zernike_res.fits';
res_phase = fitsread(res_path);
%%
figure()

plot(res_slope(:,2))
hold on;
plot(-res_phase(:,1))

title('saxo+ standalone 1st KL residual bright 1')
legend('slopes','phase')

%%
figure()
% subplot(2,2,1)
plot(f,(psd_slope(:,mode)))
hold on;
plot(f,(psd_phase(:,mode)))