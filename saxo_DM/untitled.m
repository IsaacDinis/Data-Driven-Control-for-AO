% tilt = fitsread('phase.fits');
KL = fitsread('phase2.fits');

% plop = KL/max(KL,[],'all')*max(tilt,[],'all');
% plop = KL/norm(KL,'fro')*norm(tilt,'fro');
% figure()
% imagesc(tilt+plop)

% pl = tilt+plop;
figure()
stem(KL(60,:))

figure()
imagesc(KL)