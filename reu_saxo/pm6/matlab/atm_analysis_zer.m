fs = 2760;

psd_path = '../results/bright1_03_2ms/cao/saxoplus_zernike_psd.fits';
psd_03_02ms_cao = fitsread(psd_path);

psd_path = '../results/bright1_09_2ms/cao/saxoplus_zernike_psd.fits';
psd_09_02ms_cao = fitsread(psd_path);

psd_path = '../results/bright1_09_8ms/cao/saxoplus_zernike_psd.fits';
psd_09_08ms_cao = fitsread(psd_path);

psd_path = '../results/bright1_03_8ms/cao/saxoplus_zernike_psd.fits';
psd_03_08ms_cao = fitsread(psd_path);

psd_path = '../results/bright1_03_2ms/dcao/saxoplus_zernike_psd.fits';
psd_03_02ms_dcao = fitsread(psd_path);

psd_path = '../results/bright1_09_2ms/dcao/saxoplus_zernike_psd.fits';
psd_09_02ms_dcao = fitsread(psd_path);

psd_path = '../results/bright1_09_8ms/dcao/saxoplus_zernike_psd.fits';
psd_09_08ms_dcao = fitsread(psd_path);

psd_path = '../results/bright1_03_8ms/dcao/saxoplus_zernike_psd.fits';
psd_03_08ms_dcao = fitsread(psd_path);

f = fitsread('../results/bright1_03_8ms/dcao/freq.fits');

%%
mode = 2;
figure()
% subplot(2,2,1)
semilogx(f,10*log10(psd_03_02ms_cao(:,mode)))
hold on;
semilogx(f,10*log10(psd_09_02ms_cao(:,mode)))
semilogx(f,10*log10(psd_09_08ms_cao(:,mode)))
semilogx(f,10*log10(psd_03_08ms_cao(:,mode)))
title('saxo+ standalone tilt residual PSD bright 1')
legend('0.3" 0.2ms','0.9" 0.2ms','0.9" 0.8ms','0.3" 0.8ms')
% legend_stand_g3 = sprintf('standalone: gain = 0.3 rms = %.2f',rms_standalone_g3);
% legend_dcao_g3 = sprintf('dcao: gain = 0.3 rms = %.2f',rms_dcao_g3);
% legend(legend_stand_g3,legend_dcao_g3)
% title('Tip residual, bright case, seing = 0.3", t0 = 2ms')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()


figure()
% subplot(2,2,1)
semilogx(f,10*log10(psd_03_02ms_dcao(:,mode)))
hold on;
semilogx(f,10*log10(psd_09_02ms_dcao(:,mode)))
semilogx(f,10*log10(psd_09_08ms_dcao(:,mode)))
semilogx(f,10*log10(psd_03_08ms_dcao(:,mode)))
title('saxo+ DCAO tilt residual PSD bright 1')
legend('0.3" 0.2ms','0.9" 0.2ms','0.9" 0.8ms','0.3" 0.8ms')
% legend_stand_g3 = sprintf('standalone: gain = 0.3 rms = %.2f',rms_standalone_g3);
% legend_dcao_g3 = sprintf('dcao: gain = 0.3 rms = %.2f',rms_dcao_g3);
% legend(legend_stand_g3,legend_dcao_g3)
% title('Tip residual, bright case, seing = 0.3", t0 = 2ms')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()


%%
res_path = '../results/bright1_03_2ms/cao/saxoplus_zernike_res.fits';
rms_03_02ms_cao = rms(fitsread(res_path));

res_path = '../results/bright1_03_8ms/cao/saxoplus_zernike_res.fits';
rms_03_08ms_cao = rms(fitsread(res_path));

res_path = '../results/bright1_09_2ms/cao/saxoplus_zernike_res.fits';
rms_09_02ms_cao = rms(fitsread(res_path));

res_path = '../results/bright1_09_8ms/cao/saxoplus_zernike_res.fits';
rms_09_08ms_cao = rms(fitsread(res_path));


res_path = '../results/bright1_03_2ms/dcao/saxoplus_zernike_res.fits';
rms_03_02ms_dcao = rms(fitsread(res_path));

res_path = '../results/bright1_03_8ms/dcao/saxoplus_zernike_res.fits';
rms_03_08ms_dcao = rms(fitsread(res_path));

res_path = '../results/bright1_09_2ms/dcao/saxoplus_zernike_res.fits';
rms_09_02ms_dcao = rms(fitsread(res_path));

res_path = '../results/bright1_09_8ms/dcao/saxoplus_zernike_res.fits';
rms_09_08ms_dcao = rms(fitsread(res_path));

%%

figure()
plot(rms_03_02ms_cao)
hold on;
plot(rms_09_02ms_cao)
plot(rms_09_08ms_cao)
plot(rms_03_08ms_cao)
title('saxo+ standalole zernike residual RMS bright 1 per mode')
legend('0.3" 0.2ms','0.9" 0.2ms','0.9" 0.8ms','0.3" 0.8ms')
xlabel('mode')
ylabel('RMS')

figure()
plot(rms_03_02ms_dcao)
hold on;
plot(rms_09_02ms_dcao)
plot(rms_09_08ms_dcao)
plot(rms_03_08ms_dcao)
title('saxo+ DCAO zernike residual RMS bright 1 per mode')
legend('0.3" 0.2ms','0.9" 0.2ms','0.9" 0.8ms','0.3" 0.8ms')
xlabel('mode')
ylabel('RMS')