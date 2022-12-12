DM_freq = [800,1200,1600,2000,2400,2800];

RMS_no_noise_1000 = [123,44.5,23.0,49.3,25.0,29.1,19.6]; %integrator gain = 0.5
RMS_ProxCen_1000 = [1374,1335,1320,1351,1321,1321,1284]; %integrator gain = 0.2
RMS_no_noise_1000_1delay = [17.3,12.4,11.6,12.2,11.2,10.9,15.2]; %integrator gain = 0.5

figure()
plot(DM_freq,RMS_no_noise_1000(1:end-1),'-o')
title('residual RMS, 2 frames delay, no photon noise',' RMS without DM dymanics = 19.6')
xlim([DM_freq(1),DM_freq(end)])
xlabel('DM natural frequency [Hz]')
ylabel('RMS')
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/psd.pdf -transparent

figure()
plot(DM_freq,RMS_ProxCen_1000(1:end-1),'-o')
title('residual RMS, 2 frames delay, 10.4e6 photons/m2/s ',' RMS without DM dymanics = 1284')
xlim([DM_freq(1),DM_freq(end)])
xlabel('DM natural frequency [Hz]')
ylabel('RMS')
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()

figure()
plot(DM_freq,RMS_no_noise_1000_1delay(1:end-1),'-o')
title('residual RMS, 1 frame delay, no photon noise',' RMS without DM dymanics = 15.2')
xlim([DM_freq(1),DM_freq(end)])
xlabel('DM natural frequency [Hz]')
ylabel('RMS')
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()