phase_ol = fitsread('results/dcao/ol/saxoplus_KL_phase_res.fits');
slopes_ol = fitsread('results/dcao/ol/saxoplus_KL_res.fits');
phase_cl = fitsread('results/dcao/cl/saxoplus_KL_phase_res.fits');
slopes_cl = fitsread('results/dcao/cl/saxoplus_KL_res.fits');
command_cl = fitsread('results/dcao/cl/saxoplus_KL_u.fits');

mode = 1;

figure()
plot(phase_ol(:,mode))
hold on;
plot(slopes_ol(:,mode))
plot(command_cl(:,mode))
legend('phase','slopes','command')

figure()
plot(slopes_cl(:,mode))
hold on;
plot(phase_cl(:,mode))
plot(phase_ol(:,mode)-command_cl(:,mode))

legend('slopes','phase cl','phase ol')

figure()
plot(phase_cl(:,mode)-phase_ol(:,mode)+command_cl(:,mode))