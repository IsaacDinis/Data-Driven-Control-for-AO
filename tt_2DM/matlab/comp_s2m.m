S2M_DM0 = fitsread('../compass/calib_mat/S2M_DM0.fits');
S2M_DM1 = fitsread('../compass/calib_mat/S2M_DM0.fits');
figure()
plot(S2M_DM1(1,:))
hold on
plot(S2M_DM0(1,:))

figure()
plot(S2M_DM0(2,:))
hold on
plot(S2M_DM1(2,:))