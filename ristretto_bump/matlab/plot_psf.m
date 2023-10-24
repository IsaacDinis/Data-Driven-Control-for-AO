psf_bump = fitsread('../compass/results/bump/psf.fits');
psf_no_bump = fitsread('../compass/results/no_bump/psf.fits');
figure()
imagesc(log10(psf_bump))
figure()
imagesc(log10(psf_no_bump))


stroke_bump = fitsread('../compass/results/bump/HODM_stroke.py');
stroke_no_bump = fitsread('../compass/results/no_bump/HODM_stroke.py');
figure()
plot(stroke_bump);
hold on
plot(stroke_no_bump);