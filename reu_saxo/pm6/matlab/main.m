% Bright 1
contrast_saxo = fitsread('../results/bright1_1_3ms/cao_g03/contrastCurves/saxo.fits');
contrast_saxo = contrast_saxo(2,:);

contrast_cao = fitsread('../results/bright1_1_3ms/cao_g03/contrastCurves/saxoplus.fits');
contrast_cao = contrast_cao(2,:);

contrast_dcao_int = fitsread('../results/bright1_1_3ms/dcao_gain05/contrastCurves/saxoplus.fits');
contrast_dcao_int = contrast_dcao_int(2,:);

contrast_dcao_dd = fitsread('../results/bright1_1_3ms/dcao_dd_full/contrastCurves/saxoplus.fits');
contrast_dcao_dd = contrast_dcao_dd(2,:);

figure()
semilogy(contrast_saxo)
hold on
semilogy(contrast_cao)
semilogy(contrast_dcao_int)
semilogy(contrast_dcao_dd)
legend('saxo alone','standalone','dcao integrator','dcao data-driven')
title('Radial average contrast after perfect coronograph at 1650 nm')
subtitle('Bright 1 bad seeing conditions 1.0" t0 = 3 ms')
ylabel("contrast")
xlabel("r (lamb/D)")
grid;
make_it_nicer()
set(gcf,'Position',[100 100 700 500])
%% Bright 1 good seeing
contrast_saxo = fitsread('../results/bright1_04_9ms/cao_g04/contrastCurves/saxo.fits');
contrast_saxo = contrast_saxo(2,:);

contrast_cao = fitsread('../results/bright1_04_9ms/cao_2/contrastCurves/saxoplus.fits');
contrast_cao = contrast_cao(2,:);

contrast_dcao_int = fitsread('../results/bright1_04_9ms/dcao/contrastCurves/saxoplus.fits');
contrast_dcao_int = contrast_dcao_int(2,:);

contrast_dcao_dd = fitsread('../results/bright1_04_9ms/dcao3/contrastCurves/saxoplus.fits');
contrast_dcao_dd = contrast_dcao_dd(2,:);

figure()
semilogy(contrast_saxo)
hold on
semilogy(contrast_cao)
semilogy(contrast_dcao_int)
semilogy(contrast_dcao_dd)
legend('saxo alone','standalone','dcao integrator','dcao data-driven')
title('Radial average contrast after perfect coronograph at 1650 nm')
subtitle('Bright 1 good seeing conditions 0.4" t0 = 9 ms')
ylabel("contrast")
xlabel("r (lamb/D)")
grid;
make_it_nicer()
set(gcf,'Position',[100 100 700 500])


%% Red 1
contrast_saxo = fitsread('../results/red1_04_3ms/cao/contrastCurves/saxo.fits');
contrast_saxo = contrast_saxo(2,:);

contrast_cao = fitsread('../results/red1_04_3ms/cao/contrastCurves/saxoplus.fits');
contrast_cao = contrast_cao(2,:);

contrast_dcao_int = fitsread('../results/red1_04_3ms/dcao/contrastCurves/saxoplus.fits');
contrast_dcao_int = contrast_dcao_int(2,:);

contrast_dcao_dd = fitsread('../results/red1_04_3ms/dcao_full/contrastCurves/saxoplus.fits');
contrast_dcao_dd = contrast_dcao_dd(2,:);

figure()
semilogy(contrast_saxo)
hold on
semilogy(contrast_cao)
semilogy(contrast_dcao_int)
semilogy(contrast_dcao_dd)
legend('saxo alone','standalone','dcao integrator','dcao data-driven')
title('Radial average contrast after perfect coronograph at 1650 nm')
subtitle('Red 1 seeing = 0.7" t0 = 3 ms')
ylabel("contrast")
xlabel("r (lamb/D)")
grid;
make_it_nicer()
set(gcf,'Position',[100 100 700 500])
%% Red 3 
contrast_saxo = fitsread('../results/red3_05_55ms/cao/saxo.fits');
contrast_saxo = contrast_saxo(2,:);

contrast_cao = fitsread('../results/red3_05_55ms/cao/saxoplus.fits');
contrast_cao = contrast_cao(2,:);

contrast_dcao_int = fitsread('../results/red3_05_55ms/dcao/contrastCurves/saxoplus.fits');
contrast_dcao_int = contrast_dcao_int(2,:);

contrast_dcao_dd = fitsread('../results/red3_05_55ms/dcao_dd_full/contrastCurves/saxoplus.fits');
contrast_dcao_dd = contrast_dcao_dd(2,:);

figure()
semilogy(contrast_saxo)
hold on
semilogy(contrast_cao)
semilogy(contrast_dcao_int)
semilogy(contrast_dcao_dd)
legend('saxo alone','standalone','dcao integrator','dcao data-driven')
title('Radial average contrast after perfect coronograph at 1650 nm')
subtitle('Red 3 seeing = 0.7" t0 = 5 ms')
ylabel("contrast")
xlabel("r (lamb/D)")
grid;
make_it_nicer()
set(gcf,'Position',[100 100 700 500])