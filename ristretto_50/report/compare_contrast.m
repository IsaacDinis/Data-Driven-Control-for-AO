contrast_one_dead_act = fitsread('../results/good_seeing/one_dead_act/contrast.fits');
contrast_two_dead_act_close = fitsread('../results/good_seeing/two_dead_act_close/contrast.fits');
contrast_two_dead_act_far = fitsread('../results/good_seeing/two_dead_act_far/contrast.fits');
contrast_without_defects = fitsread('../results/good_seeing/without_defects/contrast.fits');


index_LE = 2;
strehl_int_st=0;
strehl_opt_st=0;
strehl_dd_st=0;
strehl_int_dcao=0;


index_mean = 2;

avg_con_one_dead_act = mean(contrast_one_dead_act(index_mean,4:6));
avg_con_two_dead_act_close = mean(contrast_two_dead_act_close(index_mean,4:6));
avg_con_two_dead_act_far = mean(contrast_two_dead_act_far(index_mean,4:6));
avg_con_without_defects = mean(contrast_without_defects(index_mean,4:6));


figure()
index_mean = 2;
% 0.8624
semilogy(contrast_without_defects(index_mean,:))
hold on;
semilogy(contrast_one_dead_act(index_mean,:))
semilogy(contrast_two_dead_act_close(index_mean,:))
semilogy(contrast_two_dead_act_far(index_mean,:))




legend_one_dead_act = sprintf('standalone single gain int, strehl = %.2f',strehl_int_st);
legend_two_dead_act_close = sprintf('standalone modal optimized gain int, strehl = %.2f',strehl_opt_st);
legend_two_dead_act_far = sprintf('standalone order 5 datadriven, strehl = %.2f',strehl_dd_st);
legend_without_defects = sprintf('dcao single gain int, strehl = %.2f',strehl_int_dcao);



legend(legend_without_defects,legend_one_dead_act,legend_two_dead_act_close,legend_two_dead_act_far)

xlabel('r (lamb/D)')
ylabel('contrast')
title('Bright 1 best, saxo+ contrast after perfect coronograph')
subtitle("G mag = 5.5, J mag = 5.2, t0 = 9 ms, seeing = 0.4""")
grid()
set(gcf,'Position',[100 100 800 500])
make_it_nicer()