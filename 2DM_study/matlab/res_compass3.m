res_DM0_alone = fitsread('../data3/res_DM0_alone.fits')*1000;
res_tilt_DM0 = fitsread('../data3/res_tilt_DM0.fits')*1000;
res_DM0_proj = fitsread('../data3/res_DM0_proj.fits')*1000;


res_DM0_alone_4kHz = fitsread('../data3/res_DM0_alone_4kHz.fits')*1000;
res_tilt_DM0_4kHz = fitsread('../data3/res_tilt_DM0_4kHz.fits')*1000;
res_DM0_proj_4kHz = fitsread('../data3/res_DM0_proj_4kHz.fits')*1000;


t = 0:1/4000:2.5-1/4000;
t_start = 1001;
t_end = 2000;



%% RES tot DM1
% start = 30;

figure()

plot(t,res_DM0_proj)
hold on
plot(t,res_DM0_proj_4kHz)
xlabel('Time (s)')
ylabel('Res. amp. (nm)')
legend('1 kHz','4kHz','Interpreter','latex');
make_it_nicer()

set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/res_DM1.pdf -transparent -nocrop
%% RES tot DM0
figure()

plot(t,res_DM0_alone)
hold on
plot(t,res_DM0_alone_4kHz)

% plot(t,res_DM1_integrating_DM0(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp. (nm)')
title('Residual on LODM 1st KL')
legend('1 kHz','4kHz','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/res_DM0.pdf -transparent





%% RES tot phase
figure()

plot(t,res_tilt_DM0)
hold on
plot(t,res_tilt_DM0_4kHz)
% plot(t,res_DM1_integrating_DM0(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp. (um)')
title('Residual on tilt')
legend('1 kHz','4kHz','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/res_DM0.pdf -transparent


%%
n = 10;
w = 399;
fs = 4000;
[psd_tilt_DM0,f] = compute_psd(res_tilt_DM0',n,w,fs);
[psd_tilt_DM0_4kHz,f] = compute_psd(res_tilt_DM0_4kHz',n,w,fs);



figure()
loglog(f(1:end),psd_tilt_DM0(1:end))
hold on;
loglog(f(1:end),psd_tilt_DM0_4kHz(1:end))

legend('1 kHz','4kHz','Interpreter','latex');
title('Tilt residual PSD')
xlabel('Frequency (Hz)')
ylabel('Amp.')
% xlim([0,2050])
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
% export_fig ../plot/tilt_psd.pdf -transparent
%%
n = 10;
w = 399;
fs = 4000;
[psd_tilt_DM0,f] = compute_psd(res_DM0_alone',n,w,fs);
[psd_tilt_DM0_4kHz,f] = compute_psd(res_DM0_alone_4kHz',n,w,fs);



figure()
loglog(f(1:end),psd_tilt_DM0(1:end))
hold on;
loglog(f(1:end),psd_tilt_DM0_4kHz(1:end))

legend('1 kHz','4kHz','Interpreter','latex');
title('Tilt residual PSD')
xlabel('Frequency (Hz)')
ylabel('Amp.')
% xlim([0,2050])
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')