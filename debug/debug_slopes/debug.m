
slopes = fitsread('dummy/saxoplus_KL_res.fits');
phase = fitsread('dummy/saxoplus_KL_phase_res.fits');
phase_std = fitsread('dummy/saxoplus_std_phase.fits');
phase_std = phase_std.*phase_std;
phase_rms = rms(phase,2)/sqrt(121140);
slopes_rms = rms(slopes,2)/sqrt(8064);

figure()
plot(slopes_rms)
hold on;
plot(phase_rms)
plot(phase_std)

mean(phase_rms./phase_std)
mean(phase_std)