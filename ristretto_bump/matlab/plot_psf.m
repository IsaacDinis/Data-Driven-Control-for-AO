fs = 2000;
T = 2;
t = 0:1/fs:T-1/fs;
psf_bump = fitsread('../compass/results/bump_no_woof/psf.fits');
psf_no_bump = fitsread('../compass/results/no_bump_no_woof/psf.fits');
psf_bump_dyn = fitsread('../compass/results/bump_no_woof_dyn/psf.fits');
figure()
subplot(1,3,1)
imagesc(log10(psf_no_bump(450:580,450:580)))
title('no bump',' strehl = 0.757, phase rms = 63.137 nm')
subplot(1,3,2)
imagesc(log10(psf_bump_dyn(450:580,450:580)))
title('bump dyn',' strehl = 0.757, phase rms = 63.164 nm')
subplot(1,3,3)
imagesc(log10(psf_bump(450:580,450:580)))
title('static bump',' strehl = 0.754, phase rms = 63.721 nm')


stroke_bump = fitsread('../compass/results/bump_no_woof/HODM_stroke.fits');
stroke_no_bump = fitsread('../compass/results/no_bump_no_woof/HODM_stroke.fits');
stroke_bump_dyn = fitsread('../compass/results/bump_no_woof_dyn/HODM_stroke.fits');
figure()
subplot(1,2,1)
plot(t,stroke_no_bump);
hold on
plot(t,stroke_bump);
legend('no bump','static bump')
ylabel('stroke (um)')
xlabel('time (s)')

subplot(1,2,2)
plot(t,stroke_no_bump);
hold on
plot(t,stroke_bump_dyn);
legend('no bump','dyn. bump')
ylabel('stroke (um)')
xlabel('time (s)')
sgtitle('max stroke on bump neighbours actuators')









