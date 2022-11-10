res_int = load('../data/single_mode_res_int.mat').data';
res_dd = load('../data/single_mode_res_dd').data';
res_ol = load('../data/single_mode_dist').data';
% res_dd = readmatrix('/home/isaac/shesha/residuals/ghost3/residual.csv');
% res_int = readmatrix('/home/isaac/shesha/residuals/ghost3/residual_int.csv');
% rms_dd = rms(res_dd(5000:end,:));
% rms_int = rms(res_int(5000:end,:));
rms_dd = rms(res_dd);
rms_int = rms(res_int);
rms_ol = rms(res_ol);

gain = mean(100-rms_dd*100./rms_int);
figure()
plot(rms_dd)
hold on
plot(rms_int)
plot(rms_ol)
plot(rms_ol2)
title('Modal residual RMS for SPHERE in compass 1000Hz wind 8=m/s')
xlabel('mode')
ylabel('rms')
legend('Datadriven','Integrator');
% sum(res(:,2))
% sum(res(:,3))
% figure()
% subplot(2,2,1)
% plot(res(200:end,2));
% subplot(2,2,2)
% plot(res(200:end,3));
% subplot(2,2,3)
% plot(res(200:end,4));
% subplot(2,2,4)
% plot(res(200:end,5));
% legend('u','res');

fs = 1000;
window_size = 499;
n_average = 20;


[psd_dd, f] = compute_psd(res_dd,n_average,window_size,fs);
[psd_int, f] = compute_psd(res_int,n_average,window_size,fs);
[psd_ol, f] = compute_psd(res_ol,n_average,window_size,fs);

figure()
semilogx(f,10*log10(psd_dd))
hold on;
semilogx(f,10*log10(psd_int))
semilogx(f,10*log10(psd_ol))
legen('dd','int','ol')
title('atmosphere psd for tilt mode')
xlabel('frequency [Hz]')
ylabel('magnitude [dB]')
