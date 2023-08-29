phase_mode_DM0 = permute(fitsread('../data_parallel/phase_mode_DMO.fits'),[1,3,2]);
phase_mode_DM1 = permute(fitsread('../data_parallel/phase_mode_DM1.fits'),[1,3,2]);
phase_mode_tt = permute(fitsread('../data_parallel/phase_mode_tt.fits'),[1,3,2]);
phase_tt_DM = permute(fitsread('../data_parallel/phase_tt_DM.fits'),[1,3,2]);
%%
figure()
imagesc(phase_mode_DM0(:,:,2))

figure()
imagesc(phase_mode_DM1(:,:,1))

figure()
imagesc(phase_mode_tt(:,:,1))

figure()
imagesc(phase_mode_DM1(:,:,1)-phase_mode_DM0(:,:,1))
title('difference entre HODM et LODM tilt')
%%
tilt =  phase_mode_tt(:,:,1);
HODM = phase_mode_DM1(:,:,1);
norm(tilt,'fro')

fro = norm(tilt+HODM,'fro');
diff_HO_tt = fro*HODM+tilt;

figure()
imagesc(diff_HO_tt)
title('difference entre HODM et pure tilt')
colorbar

norm(diff_HO_tt,'fro')
%%
frob = norm(phase_mode_tt(:,:,1) + phase_mode_DM0(:,:,1), 'fro');
figure()
imagesc(frob*phase_mode_DM0(:,:,1)+phase_mode_tt(:,:,1))
title('difference entre LODM et pure tilt')
colorbar

%%
norm(diff_HO_tt,'fro');
LODM = phase_mode_DM0(:,:,1);
fro = norm(diff_HO_tt - LODM, 'fro');
diff = fro*0*LODM+diff_HO_tt;
norm(diff,'fro')

% figure()
% imagesc(diff)
% title('difference entre HODM et pure tilt')
% colorbar
% norm(diff)
%%
figure()
imagesc(phase_tt_DM(:,:,3))

figure()
imagesc(phase_tt_DM(:,:,1)-phase_mode_tt(:,:,1))
title('difference entre HODM et LODM tilt')

%%
tilt =  phase_mode_tt(:,:,1);
HODM = phase_tt_DM(:,:,3);
norm(HODM,'fro')
norm(tilt,'fro')
fro = norm(HODM,'fro')/norm(tilt,'fro');
norm(fro*tilt,'fro')
diff_HO_tt = HODM-fro*tilt;
norm(diff_HO_tt)
figure()
imagesc(diff_HO_tt)
title('difference entre HODM et pure tilt')
colorbar
%%
tilt =  phase_mode_tt(:,:,1);
handmade = phase_tt_DM(:,:,3);
HODM = phase_mode_DM1(:,:,1);

fro = norm(handmade,'fro')/norm(tilt,'fro');
diff_handmade_tt = handmade-fro*tilt;

fro = norm(HODM,'fro')/norm(tilt,'fro');
diff_HO_tt = HODM+fro*tilt;


figure()
subplot(1,2,1)

imagesc(diff_HO_tt)
title('difference between 1st KL and pure tilt')
subplot(1,2,2)

imagesc(diff_handmade_tt)
title('difference between handmade tilt and pure tilt')