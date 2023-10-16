fs = 2760;
dist_matrix_path = '../standalone/KL_saxoplus_res.fits';
KL_standalone_int = fitsread(dist_matrix_path);
KL_standalone_int = KL_standalone_int(:,1);
rms_KL_standalone_int = rms(KL_standalone_int);

dist_matrix_path = '../standalone_dd/KL_saxoplus_res.fits';
KL_standalone_dd = fitsread(dist_matrix_path);
KL_standalone_dd = KL_standalone_dd(:,1);
rms_KL_standalone_dd = rms(KL_standalone_dd);

dist_matrix_path = '../dcao/KL_saxoplus_res.fits';
KL_dcao_int = fitsread(dist_matrix_path);
KL_dcao_int = KL_dcao_int(:,1);
rms_KL_dcao_int = rms(KL_dcao_int);

dist_matrix_path = '../dcao_dd/KL_saxoplus_res.fits';
KL_dcao_dd = fitsread(dist_matrix_path);
KL_dcao_dd = KL_dcao_dd(:,1);
rms_KL_dcao_dd = rms(KL_dcao_dd);

dist_matrix_path = '../standalone_ol/KL_saxoplus_res.fits';
KL_standalone_ol = fitsread(dist_matrix_path);
KL_standalone_ol = KL_standalone_ol(:,1);
rms_KL_standalone_ol = rms(KL_standalone_ol);

dist_matrix_path = '../dcao_ol/KL_saxoplus_res.fits';
KL_dcao_ol = fitsread(dist_matrix_path);
KL_dcao_ol = KL_dcao_ol(:,1);
rms_KL_dcao_ol = rms(KL_dcao_ol);

dist_matrix_path = '../standalone/zernike_saxoplus_res.fits';
zer_standalone_int = fitsread(dist_matrix_path);
zer_standalone_int = zer_standalone_int(:,1);
rms_zer_standalone_int = rms(zer_standalone_int);

dist_matrix_path = '../standalone_dd/zernike_saxoplus_res.fits';
zer_standalone_dd = fitsread(dist_matrix_path);
zer_standalone_dd = zer_standalone_dd(:,1);
rms_zer_standalone_dd = rms(zer_standalone_dd);

dist_matrix_path = '../dcao/zernike_saxoplus_res.fits';
zer_dcao_int = fitsread(dist_matrix_path);
zer_dcao_int = zer_dcao_int(:,1);
rms_zer_dcao_int = rms(zer_dcao_int);

dist_matrix_path = '../dcao_dd/zernike_saxoplus_res.fits';
zer_dcao_dd = fitsread(dist_matrix_path);
zer_dcao_dd = zer_dcao_dd(:,1);
rms_zer_dcao_dd = rms(zer_dcao_dd);

dist_matrix_path = '../standalone_ol/zernike_saxoplus_res.fits';
zer_standalone_ol = fitsread(dist_matrix_path);
zer_standalone_ol = zer_standalone_ol(:,1);
rms_zer_standalone_ol = rms(zer_standalone_ol);

dist_matrix_path = '../dcao_ol/zernike_saxoplus_res.fits';
zer_dcao_ol = fitsread(dist_matrix_path);
zer_dcao_ol = zer_dcao_ol(:,1);
rms_zer_dcao_ol = rms(zer_dcao_ol);

contrast_standalone = fitsread('../standalone/contrastCurves/saxoplus.fits');
r = contrast_standalone(1,:);
contrast_standalone = contrast_standalone(2,:);

contrast_dcao = fitsread('../dcao/contrastCurves/saxoplus.fits');
contrast_dcao = contrast_dcao(2,:);
% contrast_dcao = contrast_dcao/contrast_dcao(end)*contrast_standalone(end);

%%
window_size = 100;
n_average = 20;

[psd_KL_dcao_ol, f, psd_KL_dcao_ol_std] = compute_psd(KL_dcao_ol,n_average,window_size,fs);
[psd_KL_dcao_dd, f] = compute_psd(KL_dcao_dd,n_average,window_size,fs);
[psd_KL_dcao_int, f] = compute_psd(KL_dcao_int,n_average,window_size,fs);
[psd_KL_standalone_dd, f] = compute_psd(KL_standalone_dd,n_average,window_size,fs);
[psd_KL_standalone_ol, f, psd_KL_standalone_ol_std] = compute_psd(KL_standalone_ol,n_average,window_size,fs);
[psd_KL_standalone_int, f] = compute_psd(KL_standalone_int,n_average,window_size,fs);

[psd_zer_dcao_ol, f, psd_zer_dcao_ol_std] = compute_psd(zer_dcao_ol,n_average,window_size,fs);
[psd_zer_dcao_dd, f] = compute_psd(zer_dcao_dd,n_average,window_size,fs);
[psd_zer_dcao_int, f] = compute_psd(zer_dcao_int,n_average,window_size,fs);
[psd_zer_standalone_dd, f] = compute_psd(zer_standalone_dd,n_average,window_size,fs);
[psd_zer_standalone_ol, f, psd_zer_standalone_ol_std] = compute_psd(zer_standalone_ol,n_average,window_size,fs);
[psd_zer_standalone_int, f] = compute_psd(zer_standalone_int,n_average,window_size,fs);




%%
figure()
subplot(2,2,1)
semilogx(f,10*log10(psd_KL_standalone_ol))
title('First KL standalone')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

subplot(2,2,2)
semilogx(f,10*log10(psd_KL_dcao_ol))
title('First KL dcao')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

subplot(2,2,3)
semilogx(f,10*log10(psd_zer_standalone_ol))
title('Tilt standalone')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

subplot(2,2,4)
semilogx(f,10*log10(psd_zer_dcao_ol))
title('Tilt dcao')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

sgtitle('PSD comparison open-loop')

%%

figure()
subplot(2,2,1)
errorbar(f, 10*log10(psd_KL_standalone_ol), log(10)*(psd_KL_standalone_ol_std./psd_KL_standalone_ol))
set(gca, 'XScale','log', 'YScale','linear')
title('First KL standalone')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

subplot(2,2,2)
errorbar(f, 10*log10(psd_KL_dcao_ol), log(10)*(psd_KL_dcao_ol_std./psd_KL_dcao_ol))
set(gca, 'XScale','log', 'YScale','linear')
title('First KL dcao')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

subplot(2,2,3)
errorbar(f, 10*log10(psd_zer_standalone_ol), log(10)*(psd_zer_standalone_ol_std./psd_zer_standalone_ol))
set(gca, 'XScale','log', 'YScale','linear')
title('Tilt standalone')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

subplot(2,2,4)
errorbar(f, 10*log10(psd_zer_dcao_ol), log(10)*(psd_zer_dcao_ol_std./psd_zer_dcao_ol))
set(gca, 'XScale','log', 'YScale','linear') 
title('Tilt dcao')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

sgtitle('PSD comparison open-loop')
%%
figure()
subplot(1,2,1)

semilogx(f,10*log10(psd_KL_standalone_int))
hold on
semilogx(f,10*log10(psd_KL_dcao_int))
legend(sprintf('standalone rms = %0.2f',rms_KL_standalone_int),sprintf('dcao rms = %0.2f',rms_KL_dcao_int),'Interpreter','latex')
title('First KL')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

subplot(1,2,2)
semilogx(f,10*log10(psd_zer_standalone_int))
hold on
semilogx(f,10*log10(psd_zer_dcao_int))
legend(sprintf('standalone rms = %0.2f',rms_zer_standalone_int),sprintf('dcao rms = %0.2f',rms_zer_dcao_int),'Interpreter','latex')
title('Tilt')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

sgtitle('PSD comparison')
%%
figure()

subplot(1,2,1)
semilogx(f,10*log10(psd_KL_standalone_int))
hold on
semilogx(f,10*log10(psd_KL_standalone_dd))
legend(sprintf('integrator rms = %0.2f',rms_KL_standalone_int),sprintf('data-driven rms = %0.2f',rms_KL_standalone_dd),'Interpreter','latex')
title('First KL in standalone')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

subplot(1,2,2)
semilogx(f,10*log10(psd_KL_dcao_int))
hold on
semilogx(f,10*log10(psd_KL_dcao_dd))
legend(sprintf('integrator rms = %0.2f',rms_KL_dcao_int),sprintf('data-driven rms = %0.2f',rms_KL_dcao_dd),'Interpreter','latex')
title('First KL in dcao')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

% subplot(2,2,3)
% semilogx(f,10*log10(psd_zer_standalone_int))
% hold on
% semilogx(f,10*log10(psd_zer_standalone_dd))
% legend(sprintf('integrator rms = %0.2f',rms_zer_standalone_int),sprintf('data-driven rms = %0.2f',rms_zer_standalone_dd),'Interpreter','latex')
% title('Tilt in standalone')
% xlabel('Frequency (Hz)')
% ylabel('Magnitude (dB)')
% make_it_nicer()
% 
% subplot(2,2,4)
% semilogx(f,10*log10(psd_zer_dcao_int))
% hold on
% semilogx(f,10*log10(psd_zer_dcao_dd))
% legend(sprintf('integrator rms = %0.2f',rms_zer_dcao_int),sprintf('data-driven rms = %0.2f',rms_zer_dcao_dd),'Interpreter','latex')
% title('Tilt in dcao')
% xlabel('Frequency (Hz)')
% ylabel('Magnitude (dB)')
% make_it_nicer()

sgtitle('PSD comparison')



%%

figure()

semilogy(r,contrast_standalone)
hold on
semilogy(r,contrast_dcao)
legend('standalone','dcao','Interpreter','latex')
ylabel('Contrast')
xlabel('r (lamb/d)')
title('Saxo+ contrast (perfect corono)')
make_it_nicer()

%%
% figure()
% periodogram(zer_dcao_ol)
% 
% figure()
% plot(zer_dcao_ol)
