
sine = load('../data/intersampling_sine_1000.mat').data;
ramp = load('../data/intersampling_ramp_1000.mat').data;
sine_1234 = load('../data/intersampling_sine_1234.mat').data;
ramp_1234 = load('../data/intersampling_ramp_1234.mat').data;

sine = sine(2,:);
ramp = ramp(2,:);
sine_1234 = sine_1234(2,:);
ramp_1234 = ramp_1234(2,:);

fs = 4000;
t = 0:1/fs:0.25-1/fs;
t_start = 9001;
t_end = 10000;

%%
figure()

plot(0:1/fs:0.01-1/fs,int(1:40))
% hold on
% plot(t,sine(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
title('Intersampling signal with ramp input')
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/inter_ramp.pdf -transparent

%%
figure()
subplot(1,3,1)
plot(0:1/fs:10,sine)
% hold on
% plot(t,sine(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')

make_it_nicer()

subplot(1,3,2)
plot(0:1/fs:0.01-1/fs,sine(1:40))
% hold on
% plot(t,sine(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
make_it_nicer()

subplot(1,3,3)
plot(1.56:1/fs:1.56+0.015-1/fs,sine(6251:6310))
% hold on
% plot(t,sine(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
make_it_nicer()

sgtitle('Intersampling signal with sine input (+ zoom on some time ranges)','Interpreter','latex','Fontsize',20)

set(gcf, 'Position',  [100, 100, 2100, 450])
set(gcf,'PaperType','A4')
% export_fig ../plot/inter_sine.pdf -transparent -nocrop

%%
figure()

plot(0:1/fs:0.01-1/fs,ramp_1234(1:40))
% hold on
% plot(t,ramp_1234(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
title('Intersampling signal with ramp input')
make_it_nicer()

set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')

% export_fig ../plot/inter_ramp_1234.pdf -transparent

%%
figure()

subplot(2,2,1)
plot(0:1/fs:0.025-1/fs,ramp(1:100))
xlabel('Time (s)')
ylabel('Amp.')
title('ramp input (ZOH fs = 1000 Hz)')
make_it_nicer()

subplot(2,2,2)
plot(0:1/fs:0.05-1/fs,sine(1:200))
xlabel('Time (s)')
ylabel('Amp.')
title('sine input (ZOH fs = 1000 Hz)')
make_it_nicer()

subplot(2,2,3)
plot(0:1/fs:0.025-1/fs,ramp_1234(1:100))
xlabel('Time (s)')
ylabel('Amp.')
title('ramp input (ZOH fs = 1234 Hz)')
make_it_nicer()

subplot(2,2,4)
plot(0:1/fs:0.05-1/fs,sine_1234(1:200))
xlabel('Time (s)')
ylabel('Amp.')
title('sine input (ZOH fs = 1234 Hz)')
make_it_nicer()

sgtitle('Intersampling signal','Interpreter','latex','Fontsize',20)
set(gcf, 'Position',  [100, 100, 1500, 800])
set(gcf,'PaperType','A4')

% export_fig ../plot/intersampling.pdf -transparent -nocrop

%%
figure()

plot(0:1/fs:0.05-1/fs,ramp_1234(1:200))
% hold on
% plot(t,ramp_1234(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
title('Intersampling signal with ramp input')
make_it_nicer()

set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')

% export_fig ../plot/inter_ramp_1234.pdf -transparent
%% PSD
n = 80;
w = 500;

[inter_psd_sine,f] = compute_psd((sine-mean(sine))',n,w,fs);
[inter_psd_ramp,f] = compute_psd((int-mean(int))',n,w,fs);
[inter_psd_sine_1234,f] = compute_psd((sine_1234-mean(sine_1234))',n,w,fs);
[inter_psd_ramp_1234,f] = compute_psd((ramp_1234-mean(ramp_1234))',n,w,fs);

figure()
subplot(2,2,1)
semilogy(f(1:end),inter_psd_ramp(1:end))
title('ramp input (ZOH fs = 1000 Hz)')
xlabel('Frequency (Hz)')
ylabel('Amp.')
xlim([0,2050])
make_it_nicer()

subplot(2,2,2)
semilogy(f(1:end),inter_psd_sine(1:end))
title('sine input (ZOH fs = 1000 Hz)')
xlabel('Frequency (Hz)')
ylabel('Amp.')
xlim([0,2050])
make_it_nicer()

subplot(2,2,3)
semilogy(f(1:end),inter_psd_ramp_1234(1:end))
title('ramp input (ZOH fs = 1234 Hz)')
xlabel('Frequency (Hz)')
ylabel('Amp.')
xlim([0,2050])
make_it_nicer()

subplot(2,2,4)
semilogy(f(1:end),inter_psd_sine_1234(1:end))
title('sine input (ZOH fs = 1234 Hz)')
xlabel('Frequency (Hz)')
ylabel('Amp.')
xlim([0,2050])
make_it_nicer()
sgtitle('Intersampling PSD','Interpreter','latex','Fontsize',20)
set(gcf, 'Position',  [100, 100, 1500, 800])
set(gcf,'PaperType','A4')

% export_fig ../plot/intersampling_psd.pdf -transparent -nocrop