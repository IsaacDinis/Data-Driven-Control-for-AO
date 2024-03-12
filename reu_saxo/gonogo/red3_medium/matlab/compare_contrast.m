contrast_int = fitsread('../results/standalone/integrator_02/contrastCurves/saxoplus.fits');
contrast_opt = fitsread('../results/standalone/modal_gain_opt/contrastCurves/saxoplus.fits');
contrast_dd = fitsread('../results/standalone/datadriven_5/contrastCurves/saxoplus.fits');

strehl_int = fitsread('../results/standalone/integrator_02/strehlCurves/saxoplus.fits');
strehl_opt = fitsread('../results/standalone/modal_gain_opt/strehlCurves/saxoplus.fits');
strehl_dd = fitsread('../results/standalone/datadriven_5/strehlCurves/saxoplus.fits');

index_LE = 2;
strehl_int=strehl_int(index_LE,end)*100;
strehl_opt=strehl_opt(index_LE,end)*100;
strehl_dd=strehl_dd(index_LE,end)*100;

figure()
index_mean = 2;
% 0.8624
semilogy(contrast_int(index_mean,:))
hold on;
semilogy(contrast_opt(index_mean,:))
semilogy(contrast_dd(index_mean,:))

legend_int = sprintf('single gain int, strehl = %.2f',strehl_int);
legend_opt = sprintf('modal optimized gain int, strehl = %.2f',strehl_opt);
legend_dd = sprintf('order 5 datadriven, strehl = %.2f',strehl_dd);

legend(legend_int,legend_opt,legend_dd)

xlabel('r (lamb/D)')
ylabel('contrast')
title('Red 3 medium, saxo+ contrast after perfect coronograph')
subtitle("G mag = 14.5, J mag = 10.1, t0 = 5.5 ms, seeing = 0.7""")
grid()
make_it_nicer()