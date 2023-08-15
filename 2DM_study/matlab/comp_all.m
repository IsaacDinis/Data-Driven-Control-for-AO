res_DM0_all = fitsread('../data_parallel/res_DM0_all.fits');
res_DM1_all = fitsread('../data_parallel/res_DM1_all.fits');
res_phase = fitsread('../data_parallel/res_phase.fits');

res_DM0_all_DM1_alone = fitsread('../data_parallel/res_DM0_all_DM1_alone.fits');
res_DM1_all_DM1_alone = fitsread('../data_parallel/res_DM1_all_DM1_alone.fits');
res_phase_DM1_alone = fitsread('../data_parallel/res_phase_DM1_alone.fits');

res_DM0_all_openloop = fitsread('../data_parallel/res_DM0_all_openloop.fits');
res_DM1_all_openloop = fitsread('../data_parallel/res_DM1_all_openloop.fits');
res_phase_openloop = fitsread('../data_parallel/res_phase_openloop.fits');

voltage_DM1_all = fitsread('../data_parallel/voltage_DM1_all.fits');
voltage_DM1_all_DM1_alone = fitsread('../data_parallel/voltage_DM1_all_DM1_alone.fits');

fs = 1000;

t = 0:1/fs:4.99-1/fs;
%%
res_DM0_all_rms = rms(res_DM0_all(:,150:end),2);
res_DM1_all_rms = rms(res_DM1_all(:,150:end),2);
res_phase_rms = rms(res_phase(:,150:end),2);

res_DM0_all_DM1_alone_rms = rms(res_DM0_all_DM1_alone(:,150:end),2);
res_DM1_all_DM1_alone_rms = rms(res_DM1_all_DM1_alone(:,150:end),2);
res_phase_DM1_alone_rms = rms(res_phase_DM1_alone(:,150:end),2);

res_phase_openloop_rms = rms(res_phase_openloop(:,150:end),2);

voltage_DM1_all_rms = rms(voltage_DM1_all(:,150:end),2);
voltage_DM1_all_DM1_alone_rms = rms(voltage_DM1_all_DM1_alone(:,150:end),2);

%%
figure()
stem(res_DM0_all_rms)
hold on
stem(res_DM0_all_DM1_alone_rms)
title('LODM')
xlabel('mode')
legend('both DMs','HODM only')
ylabel('rms')

figure()
stem(res_DM1_all_rms)
hold on
stem(res_DM1_all_DM1_alone_rms)
title('HODM')
xlabel('mode')
legend('both DMs','HODM only')
ylabel('rms')

figure()
stem(res_phase_rms)
hold on
stem(res_phase_DM1_alone_rms)
title('phase')
xlabel('mode')
legend('both DMs','HODM only')
ylabel('rms')

figure()
stem(res_phase_openloop_rms)
title('phase openloop')
xlabel('mode')
ylabel('rms')

figure()
stem(voltage_DM1_all_rms)
hold on
stem(voltage_DM1_all_DM1_alone_rms)
title('DM1 voltage')
legend('both DMs','HODM only')
xlabel('mode')
ylabel('rms')