fs = 2760;
dist_matrix_path = '../standalone/KL_saxoplus_res.fits';
KL_standalone_int = fitsread(dist_matrix_path);
KL_standalone_int = KL_standalone_int(:,1);

dist_matrix_path = '../standalone_dd/KL_saxoplus_res.fits';
KL_standalone_dd = fitsread(dist_matrix_path);
KL_standalone_dd = KL_standalone_dd(:,1);

dist_matrix_path = '../dcao/KL_saxoplus_res.fits';
KL_dcao_int = fitsread(dist_matrix_path);
KL_dcao_int = KL_dcao_int(:,1);

dist_matrix_path = '../dcao_dd/KL_saxoplus_res.fits';
KL_dcao_dd = fitsread(dist_matrix_path);
KL_dcao_dd = KL_dcao_dd(:,1);

dist_matrix_path = '../standalone_ol/KL_saxoplus_res.fits';
KL_standalone_ol = fitsread(dist_matrix_path);
KL_standalone_ol = KL_standalone_ol(:,1);

dist_matrix_path = '../dcao_ol/KL_saxoplus_res.fits';
KL_dcao_ol = fitsread(dist_matrix_path);
KL_dcao_ol = KL_dcao_ol(:,1);

dist_matrix_path = '../standalone/zernike_saxoplus_res.fits';
zer_standalone_int = fitsread(dist_matrix_path);
zer_standalone_int = zer_standalone_int(:,1);

dist_matrix_path = '../standalone_dd/zernike_saxoplus_res.fits';
zer_standalone_dd = fitsread(dist_matrix_path);
zer_standalone_dd = zer_standalone_dd(:,1);

dist_matrix_path = '../dcao/zernike_saxoplus_res.fits';
zer_dcao_int = fitsread(dist_matrix_path);
zer_dcao_int = zer_dcao_int(:,1);

dist_matrix_path = '../dcao_dd/zernike_saxoplus_res.fits';
zer_dcao_dd = fitsread(dist_matrix_path);
zer_dcao_dd = zer_dcao_dd(:,1);

dist_matrix_path = '../standalone_ol/zernike_saxoplus_res.fits';
zer_standalone_ol = fitsread(dist_matrix_path);
zer_standalone_ol = zer_standalone_ol(:,1);

dist_matrix_path = '../dcao_ol/zernike_saxoplus_res.fits';
zer_dcao_ol = fitsread(dist_matrix_path);
zer_dcao_ol = zer_dcao_ol(:,1);

contrast_standalone = fitsread('../standalone/contrastCurves/saxoplus.fits');
r = contrast_standalone(1,:);
contrast_standalone = contrast_standalone(2,:);

contrast_dcao = fitsread('../dcao/contrastCurves/saxoplus.fits');
contrast_dcao = contrast_dcao(2,:);
% contrast_dcao = contrast_dcao/contrast_dcao(end)*contrast_standalone(end);

%%
window_size = 200;
n_average = 10;

[psd_KL_dcao_ol, f] = compute_psd(KL_dcao_ol,n_average,window_size,fs);
[psd_KL_dcao_dd, f] = compute_psd(KL_dcao_dd,n_average,window_size,fs);
[psd_KL_dcao_int, f] = compute_psd(KL_dcao_int,n_average,window_size,fs);
[psd_KL_standalone_dd, f] = compute_psd(KL_standalone_dd,n_average,window_size,fs);
[psd_KL_standalone_ol, f] = compute_psd(KL_standalone_ol,n_average,window_size,fs);
[psd_KL_standalone_int, f] = compute_psd(KL_standalone_int,n_average,window_size,fs);

[psd_zer_dcao_ol, f] = compute_psd(zer_dcao_ol,n_average,window_size,fs);
[psd_zer_dcao_dd, f] = compute_psd(zer_dcao_dd,n_average,window_size,fs);
[psd_zer_dcao_int, f] = compute_psd(zer_dcao_int,n_average,window_size,fs);
[psd_zer_standalone_dd, f] = compute_psd(zer_standalone_dd,n_average,window_size,fs);
[psd_zer_standalone_ol, f] = compute_psd(zer_standalone_ol,n_average,window_size,fs);
[psd_zer_standalone_int, f] = compute_psd(zer_standalone_int,n_average,window_size,fs);




%%
figure()
subplot(2,2,1)
semilogx(f,10*log10(psd_KL_standalone_ol))
title('First KL standalone')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

subplot(2,2,2)
semilogx(f,10*log10(psd_KL_dcao_ol))
title('First KL dcao')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

subplot(2,2,3)
semilogx(f,10*log10(psd_zer_standalone_ol))
title('Tilt standalone')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

subplot(2,2,4)
semilogx(f,10*log10(psd_zer_dcao_ol))
title('Tilt dcao')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

sgtitle('PSD comparison open-loop')
%%
figure()
subplot(1,2,1)

semilogx(f,10*log10(psd_KL_standalone_int))
hold on
semilogx(f,10*log10(psd_KL_dcao_int))
legend('standalone','dcao')
title('First KL')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

subplot(1,2,2)
semilogx(f,10*log10(psd_zer_standalone_int))
hold on
semilogx(f,10*log10(psd_zer_dcao_int))
legend('standalone','dcao')
title('Tilt')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

sgtitle('PSD comparison')
%%
figure()

subplot(2,2,1)
semilogx(f,10*log10(psd_KL_standalone_int))
hold on
semilogx(f,10*log10(psd_KL_standalone_dd))
legend('integrator','data-driven')
title('First KL in standalone')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

subplot(2,2,2)
semilogx(f,10*log10(psd_KL_dcao_int))
hold on
semilogx(f,10*log10(psd_KL_dcao_dd))
legend('integrator','data-driven')
title('First KL in dcao')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

subplot(2,2,3)
semilogx(f,10*log10(psd_zer_standalone_int))
hold on
semilogx(f,10*log10(psd_zer_standalone_dd))
legend('integrator','data-driven')
title('Tilt in standalone')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

subplot(2,2,4)
semilogx(f,10*log10(psd_zer_dcao_int))
hold on
semilogx(f,10*log10(psd_zer_dcao_dd))
legend('integrator','data-driven')
title('Tilt in dcao')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

sgtitle('PSD comparison')



%%

figure()

semilogy(r,contrast_standalone)
hold on
semilogy(r,contrast_dcao)
legend('standalone','dcao')
ylabel('Contrast')
xlabel('r (lamb/d)')
title('Saxo+ contrast (perfect corono)')




