res_DM0_all = fitsread('../data_parallel/res_DM0_all.fits');
res_DM1_all = fitsread('../data_parallel/res_DM1_all.fits');
res_phase = fitsread('../data_parallel/res_phase.fits');

res_DM0_all_DM1_alone = fitsread('../data_parallel/res_DM0_all_DM1_alone.fits');
res_DM1_all_DM1_alone = fitsread('../data_parallel/res_DM1_all_DM1_alone.fits');
res_phase_DM1_alone = fitsread('../data_parallel/res_phase_DM1_alone.fits');

res_DM0_all_openloop = fitsread('../data_parallel/res_DM0_all_openloop.fits');
res_DM1_all_openloop = fitsread('../data_parallel/res_DM1_all_openloop.fits');
res_phase_openloop = fitsread('../data_parallel/res_phase_openloop.fits');

fs = 1000;

t = 0:1/fs:4.99-1/fs;
%%
res_DM0_all_rms = rms(res_DM0_all,2);
res_DM1_all_rms = rms(res_DM1_all,2);
res_phase_rms = rms(res_phase,2);

res_DM0_all_DM1_alone_rms = rms(res_DM0_all_DM1_alone,2);
res_DM1_all_DM1_alone_rms = rms(res_DM1_all_DM1_alone,2);
res_phase_DM1_alone_rms = rms(res_phase_DM1_alone,2);

res_phase_openloop_rms = rms(res_phase_openloop,2);
%%
figure()
stem(res_DM0_all_rms)
hold on
stem(res_DM0_all_DM1_alone_rms)
title('LODM')
legend('both DMs','HODM only')

figure()
stem(res_DM1_all_rms)
hold on
stem(res_DM1_all_DM1_alone_rms)
title('HODM')
legend('both DMs','HODM only')

figure()
stem(res_phase_rms)
hold on
stem(res_phase_DM1_alone_rms)
title('phase')
legend('both DMs','HODM only')

figure()
stem(res_phase_openloop_rms)
title('phase openloop')

