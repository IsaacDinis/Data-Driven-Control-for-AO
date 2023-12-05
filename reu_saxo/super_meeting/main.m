fs = 2760;
mode = 3;

dist_matrix_path = 'results/standalone/bright1_s0_3_t2_g3/saxoplus_KL_res_psd.fits';
KL_standalone_psd_g3 = fitsread(dist_matrix_path);
f = KL_standalone_psd_g3(:,1);
KL_standalone_psd_g3 = KL_standalone_psd_g3(:,mode);
rms_standalone_g3 = sqrt(sum(KL_standalone_psd_g3)*f(2)*275/fs);

dist_matrix_path = 'results/standalone/bright1_s0_3_t2_g2/saxoplus_KL_res_psd.fits';
KL_standalone_psd_g2 = fitsread(dist_matrix_path);
KL_standalone_psd_g2 = KL_standalone_psd_g2(:,mode);
rms_standalone_g2 = sqrt(sum(KL_standalone_psd_g2)*f(2)*275/fs);

dist_matrix_path = 'results/dcao/bright1_s0_3_t2_g3/saxoplus_KL_res_psd.fits';
KL_dcao_psd_g3 = fitsread(dist_matrix_path);
KL_dcao_psd_g3 = KL_dcao_psd_g3(:,mode);
rms_dcao_g3 = sqrt(sum(KL_dcao_psd_g3)*f(2)*275/fs);

dist_matrix_path = 'results/dcao/bright1_s0_3_t2_g5/saxoplus_KL_res_psd.fits';
KL_dcao_psd_g5 = fitsread(dist_matrix_path);
KL_dcao_psd_g5 = KL_dcao_psd_g5(:,mode);
rms_dcao_g5 = sqrt(sum(KL_dcao_psd_g5)*f(2)*275/fs);


%%
figure()
% subplot(2,2,1)
semilogx(f,10*log10(KL_standalone_psd_g3))
hold on;
semilogx(f,10*log10(KL_dcao_psd_g3))
legend_stand_g3 = sprintf('standalone: gain = 0.3 rms = %.2f',rms_standalone_g3);
legend_dcao_g3 = sprintf('dcao: gain = 0.3 rms = %.2f',rms_dcao_g3);
legend(legend_stand_g3,legend_dcao_g3)
title('Tip residual, bright case, seing = 0.3", t0 = 2ms')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

figure()
% subplot(2,2,1)
semilogx(f,10*log10(KL_dcao_psd_g3))
hold on;
semilogx(f,10*log10(KL_dcao_psd_g5))
legend_dcao_g3 = sprintf('dcao: gain = 0.3 rms = %.2f',rms_dcao_g3);
legend_dcao_g5 = sprintf('dcao: gain = 0.5 rms = %.2f',rms_dcao_g5);
legend(legend_dcao_g3,legend_dcao_g5)
title('Tip residual, bright case, seing = 0.3", t0 = 2ms')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

figure()
% subplot(2,2,1)
semilogx(f,10*log10(KL_standalone_psd_g3))
hold on;
semilogx(f,10*log10(KL_standalone_psd_g2))
legend_stand_g3 = sprintf('standalone: gain = 0.3 rms = %.2f',rms_standalone_g3);
legend_stand_g2 = sprintf('standalone: gain = 0.2 rms = %.2f',rms_standalone_g2);
legend(legend_stand_g3,legend_stand_g2)
title('Tip residual, bright case, seing = 0.3", t0 = 2ms')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

