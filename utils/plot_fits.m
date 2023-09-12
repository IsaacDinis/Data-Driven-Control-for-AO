
info = fitsinfo('HODM_gauss_fitSPARTA_deadactus.fits');
disp(info.Contents);

DM1 = fitsread('HODM_gauss_fitSPARTA_deadactus.fits','Image',1);
DM2 = fitsread('HODM_gauss_fitSPARTA_deadactus.fits','Image',2);
DM3 = fitsread('HODM_gauss_fitSPARTA_deadactus.fits','Image',3);

figure()
imagesc(squeeze(DM2(:,1,:)))

plop = max(max(DM2,[],1),[],3);
dummy = find(plop ==0);
figure()
stem(plop)