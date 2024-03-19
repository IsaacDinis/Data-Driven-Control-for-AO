contrast_int_st = fitsread('../results/standalone/integrator_02/contrastCurves/saxoplus.fits');
contrast_opt_st = fitsread('../results/standalone/modal_gain_opt/contrastCurves/saxoplus.fits');
contrast_dd_st = fitsread('../results/standalone/datadriven_5/contrastCurves/saxoplus.fits');

contrast_int_dcao = fitsread('../results/dcao/integrator/contrastCurves/saxoplus.fits');
contrast_opt_dcao = fitsread('../results/dcao/modal_gain_opt/contrastCurves/saxoplus.fits');
contrast_dd_dcao = fitsread('../results/dcao/datadriven_5/contrastCurves/saxoplus.fits');

strehl_int_st = fitsread('../results/standalone/integrator_02/strehlCurves/saxoplus.fits');
strehl_opt_st = fitsread('../results/standalone/modal_gain_opt/strehlCurves/saxoplus.fits');
strehl_dd_st = fitsread('../results/standalone/datadriven_5/strehlCurves/saxoplus.fits');

strehl_int_dcao = fitsread('../results/dcao/integrator/strehlCurves/saxoplus.fits');
strehl_opt_dcao = fitsread('../results/dcao/modal_gain_opt/strehlCurves/saxoplus.fits');
strehl_dd_dcao = fitsread('../results/dcao/datadriven_5/strehlCurves/saxoplus.fits');

index_LE = 2;
strehl_int_st=strehl_int_st(index_LE,end)*100;
strehl_opt_st=strehl_opt_st(index_LE,end)*100;
strehl_dd_st=strehl_dd_st(index_LE,end)*100;

strehl_int_dcao=strehl_int_dcao(index_LE,end)*100;
strehl_opt_dcao=strehl_opt_dcao(index_LE,end)*100;
strehl_dd_dcao=strehl_dd_dcao(index_LE,end)*100;

index_mean = 2;

avg_con_int_st = mean(contrast_int_st(index_mean,4:6));
avg_con_opt_st = mean(contrast_opt_st(index_mean,4:6));
avg_con_dd_st = mean(contrast_dd_st(index_mean,4:6));

avg_con_int_dcao = mean(contrast_int_dcao(index_mean,4:6));
avg_con_opt_dcao = mean(contrast_opt_dcao(index_mean,4:6));
avg_con_dd_dcao = mean(contrast_dd_dcao(index_mean,4:6));

figure()
index_mean = 2;
% 0.8624
semilogy(contrast_int_st(index_mean,:),'c')
hold on;
semilogy(contrast_opt_st(index_mean,:),'r')
semilogy(contrast_dd_st(index_mean,:),'b')
semilogy(contrast_int_dcao(index_mean,:),'c--')
semilogy(contrast_opt_dcao(index_mean,:),'r--')
semilogy(contrast_dd_dcao(index_mean,:),'b--')

legend_int_st = sprintf('standalone single gain int, strehl = %.2f',strehl_int_st);
legend_opt_st = sprintf('standalone modal optimized gain int, strehl = %.2f',strehl_opt_st);
legend_dd_st = sprintf('standalone order 5 datadriven, strehl = %.2f',strehl_dd_st);
legend_int_dcao = sprintf('dcao single gain int, strehl = %.2f',strehl_int_dcao);
legend_opt_dcao = sprintf('dcao modal optimized gain int, strehl = %.2f',strehl_opt_dcao);
legend_dd_dcao = sprintf('dcao order 5 datadriven, strehl = %.2f',strehl_dd_dcao);


legend(legend_int_st,legend_opt_st,legend_dd_st,legend_int_dcao,legend_opt_dcao,legend_dd_dcao)

xlabel('r (lamb/D)')
ylabel('contrast')
title('Bright 1 best, saxo+ contrast after perfect coronograph')
subtitle("G mag = 5.5, J mag = 5.2, t0 = 9 ms, seeing = 0.4""")
grid()
set(gcf,'Position',[100 100 800 500])
make_it_nicer()