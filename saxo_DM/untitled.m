tilt = fitsread('tilt.fits');
KL = fitsread('KL0.fits');

plop = KL/max(KL,[],'all')*max(tilt,[],'all');
plop = KL/norm(KL,'fro')*norm(tilt,'fro');
figure()
imagesc(tilt+plop)

pl = tilt+plop;
figure()
plot(KL(:,60))

figure()
imagesc(plop)