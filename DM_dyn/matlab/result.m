% DM_freq = [800,1200,1600,2000,2400,2800];
% 
% RMS_no_noise_1000 = [123,44.5,23.0,49.3,25.0,29.1,19.6]; %integrator gain = 0.5
% RMS_ProxCen_1000 = [1374,1335,1320,1351,1321,1321,1284]; %integrator gain = 0.2
% RMS_no_noise_1000_1delay = [17.3,12.4,11.6,12.2,11.2,10.9,15.2]; %integrator gain = 0.5
% 
% figure()
% plot(DM_freq,RMS_no_noise_1000(1:end-1),'-o')
% title('residual RMS, 2 frames delay, no photon noise',' RMS without DM dymanics = 19.6')
% xlim([DM_freq(1),DM_freq(end)])
% xlabel('DM natural frequency [Hz]')
% ylabel('RMS')
% make_it_nicer()
% set(gcf, 'Position',  [100, 100, 700, 450])
% set(gcf,'PaperType','A4')
% make_it_nicer()
% % export_fig ../plot/psd.pdf -transparent
% 
% figure()
% plot(DM_freq,RMS_ProxCen_1000(1:end-1),'-o')
% title('residual RMS, 2 frames delay, 10.4e6 photons/m2/s ',' RMS without DM dymanics = 1284')
% xlim([DM_freq(1),DM_freq(end)])
% xlabel('DM natural frequency [Hz]')
% ylabel('RMS')
% make_it_nicer()
% set(gcf, 'Position',  [100, 100, 700, 450])
% set(gcf,'PaperType','A4')
% make_it_nicer()
% 
% figure()
% plot(DM_freq,RMS_no_noise_1000_1delay(1:end-1),'-o')
% title('residual RMS, 1 frame delay, no photon noise',' RMS without DM dymanics = 15.2')
% xlim([DM_freq(1),DM_freq(end)])
% xlabel('DM natural frequency [Hz]')
% ylabel('RMS')
% make_it_nicer()
% set(gcf, 'Position',  [100, 100, 700, 450])
% set(gcf,'PaperType','A4')
% make_it_nicer()

%%
DM_freq = [700,800,850,900,950,1000,1200,1400,1600,1800,2000,2200,2400,2800,...
    3000,3400,4000];

load('out_2frame.mat')
RMS_no_noise_1000 = [out.h700,out.h800,inf,inf,out.h950,out.h1000...
    ,out.h1200,out.h1400,out.h1600,out.h1800,out.h2000,out.h2200,out.h2400...
    ,out.h2800,out.h3000,out.h3400,out.h4000,out.no_dm]; %integrator gain = 0.5

load('out_noise.mat')
RMS_ProxCen_1000 = [out.h700,out.h800,out.h850,out.h900,out.h950,out.h1000...
    ,out.h1200,out.h1400,out.h1600,out.h1800,out.h2000,out.h2200,out.h2400...
    ,out.h2800,out.h3000,out.h3400,out.h4000,out.no_dm]; %integrator gain = 0.5

load('out_1frame.mat')
RMS_no_noise_1000_1delay = [out.h700,out.h800,out.h850,out.h900,out.h950,out.h1000...
    ,out.h1200,out.h1400,out.h1600,out.h1800,out.h2000,out.h2200,out.h2400...
    ,out.h2800,out.h3000,out.h3400,out.h4000,out.no_dm]; %integrator gain = 0.5

figure()
plot(DM_freq,RMS_no_noise_1000(1:end-1),'-o')
hold on;
yl = yline(RMS_no_noise_1000(end),'--r','RMS w/o DM dyn.');
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'center';

title('Residual RMS, 2 frames delay, no photon noise')
xlim([DM_freq(1),DM_freq(end)])
xlabel('DM natural frequency (Hz)')
ylabel('RMS')
ylim([0,60]);
grid on
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/result_2frame.pdf -transparent

figure()
plot(DM_freq,RMS_ProxCen_1000(1:end-1),'-o')
hold on;
yl = yline(RMS_ProxCen_1000(end),'--r','RMS w/o DM dyn.');
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'center';
title('Residual RMS, 2 frames delay, 10.4e6 photons/m2/s ')
xlim([DM_freq(1),DM_freq(end)])
xlabel('DM natural frequency (Hz)')
ylabel('RMS')
make_it_nicer()
grid on
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/result_noise.pdf -transparent

figure()
plot(DM_freq,RMS_no_noise_1000_1delay(1:end-1),'-o')
hold on;
yl = yline(RMS_no_noise_1000_1delay(end),'--r','RMS w/o DM dyn.');
yl.LabelVerticalAlignment = 'middle';
yl.LabelHorizontalAlignment = 'center';
title('Residual RMS, 1 frame delay, no photon noise')
xlim([DM_freq(1),DM_freq(end)])
xlabel('DM natural frequency (Hz)')
ylabel('RMS')
make_it_nicer()
grid on
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')

% export_fig ../plot/result_1frame.pdf -transparent