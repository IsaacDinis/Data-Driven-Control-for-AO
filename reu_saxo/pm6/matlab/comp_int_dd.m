fs = 2760;

psd_path = '../results/bright1_03_8ms/cao/saxoplus_KL_psd.fits';
psd_int_cao = fitsread(psd_path);

psd_path = '../results/bright1_03_8ms/cao_dd/saxoplus_KL_psd.fits';
psd_dd_cao = fitsread(psd_path);

psd_path = '../results/bright1_03_8ms/dcao/saxoplus_KL_psd.fits';
psd_int_dcao = fitsread(psd_path);

psd_path = '../results/bright1_03_8ms/dcao_dd/saxoplus_KL_psd.fits';
psd_dd_dcao = fitsread(psd_path);

%%
res_path = '../results/bright1_03_8ms/cao/saxoplus_KL_res.fits';
res_int_cao = fitsread(res_path);

res_path = '../results/bright1_03_8ms/cao_dd/saxoplus_KL_res.fits';
res_dd_cao = fitsread(res_path);

res_path = '../results/bright1_03_8ms/dcao/saxoplus_KL_res.fits';
res_int_dcao = fitsread(res_path);

res_path = '../results/bright1_03_8ms/dcao_dd/saxoplus_KL_res.fits';
res_dd_dcao = fitsread(res_path);

f = fitsread('../results/bright1_03_8ms/dcao/freq.fits');

%%

mode = 1;

rms_dd_cao = rms(res_dd_cao(:,mode));
rms_int_cao = rms(res_int_cao(:,mode));
rms_dd_dcao = rms(res_dd_dcao(:,mode));
rms_int_dcao = rms(res_int_dcao(:,mode));

figure()
% subplot(2,2,1)
semilogx(f,10*log10(psd_int_cao(:,mode)))
hold on;
semilogx(f,10*log10(psd_dd_cao(:,mode)))
semilogx(f,10*log10(psd_int_dcao(:,mode)))
semilogx(f,10*log10(psd_dd_dcao(:,mode)))

title('saxo+ 1st KL residual PSD bright 1 seing 0.3" t0 8ms')
legend_dd_dcao = sprintf('dcao dd rms = %.2f',rms_dd_dcao);
legend_int_dcao = sprintf('dcao int rms = %.2f',rms_int_dcao);
legend_dd_cao = sprintf('cao dd rms = %.2f',rms_dd_cao);
legend_int_cao = sprintf('cao int rms = %.2f',rms_int_cao);


legend(legend_int_cao,legend_dd_cao,legend_int_dcao,legend_dd_dcao)
% legend_stand_g3 = sprintf('standalone: gain = 0.3 rms = %.2f',rms_standalone_g3);
% legend_dcao_g3 = sprintf('dcao: gain = 0.3 rms = %.2f',rms_dcao_g3);
% legend(legend_stand_g3,legend_dcao_g3)
% title('Tip residual, bright case, seing = 0.3", t0 = 2ms')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

%%
psd_path = '../results/bright1_03_8ms/cao/saxoplus_KL_phase_psd.fits';
psd_int_cao = fitsread(psd_path);

psd_path = '../results/bright1_03_8ms/cao_dd/saxoplus_KL_phase_psd.fits';
psd_dd_cao = fitsread(psd_path);

psd_path = '../results/bright1_03_8ms/dcao/saxoplus_KL_phase_psd.fits';
psd_int_dcao = fitsread(psd_path);

psd_path = '../results/bright1_03_8ms/dcao_dd/saxoplus_KL_phase_psd.fits';
psd_dd_dcao = fitsread(psd_path);

%%
res_path = '../results/bright1_03_8ms/cao/saxoplus_KL_phase_res.fits';
res_int_cao = fitsread(res_path);

res_path = '../results/bright1_03_8ms/cao_dd/saxoplus_KL_phase_res.fits';
res_dd_cao = fitsread(res_path);

res_path = '../results/bright1_03_8ms/dcao/saxoplus_KL_phase_res.fits';
res_int_dcao = fitsread(res_path);

res_path = '../results/bright1_03_8ms/dcao_dd/saxoplus_KL_phase_res.fits';
res_dd_dcao = fitsread(res_path);

f = fitsread('../results/bright1_03_8ms/dcao/freq.fits');

%%

mode = 1;

rms_dd_cao = rms(res_dd_cao(:,mode));
rms_int_cao = rms(res_int_cao(:,mode));
rms_dd_dcao = rms(res_dd_dcao(:,mode));
rms_int_dcao = rms(res_int_dcao(:,mode));

figure()
% subplot(2,2,1)
semilogx(f,10*log10(psd_int_cao(:,mode)))
hold on;
semilogx(f,10*log10(psd_dd_cao(:,mode)))
semilogx(f,10*log10(psd_int_dcao(:,mode)))
semilogx(f,10*log10(psd_dd_dcao(:,mode)))

title('saxo+ 1st KL residual PSD bright 1 seing 0.3" t0 8ms')
legend_dd_dcao = sprintf('dcao dd rms = %.2f',rms_dd_dcao);
legend_int_dcao = sprintf('dcao int rms = %.2f',rms_int_dcao);
legend_dd_cao = sprintf('cao dd rms = %.2f',rms_dd_cao);
legend_int_cao = sprintf('cao int rms = %.2f',rms_int_cao);


legend(legend_int_cao,legend_dd_cao,legend_int_dcao,legend_dd_dcao)
% legend_stand_g3 = sprintf('standalone: gain = 0.3 rms = %.2f',rms_standalone_g3);
% legend_dcao_g3 = sprintf('dcao: gain = 0.3 rms = %.2f',rms_dcao_g3);
% legend(legend_stand_g3,legend_dcao_g3)
% title('Tip residual, bright case, seing = 0.3", t0 = 2ms')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

%%
% figure()
% 
% plot(res_int_cao(:,mode))
% hold on;
% plot(res_dd_cao(:,mode))
% 
% title('saxo+ standalone 1st KL residual bright 1')
% legend('int','dd')
% rms(res_int_cao(:,mode))
% rms(res_dd_cao(:,mode))
%%
res_path = '../results/bright1_1_3ms/cao/saxoplus_KL_res.fits';
res_int_cao = fitsread(res_path);

res_path = '../results/bright1_04_9ms/dcao/saxoplus_KL_res.fits';
res_int_dcao = fitsread(res_path);

res_path = '../results/bright1_04_9ms/dcao_dd/saxoplus_KL_res.fits';
res_dd_dcao = fitsread(res_path);

mode = 1;


rms_int_cao = rms(res_int_cao(1:end,mode));
rms_dd_dcao = rms(res_dd_dcao(1:end,mode));
rms_int_dcao = rms(res_int_dcao(1:end,mode));

figure()
plot(res_int_cao(1:1000,mode));
hold on;
plot(res_int_dcao(1:1000,mode));
plot(res_dd_dcao(1:1000,mode));