contrast_no_dyn = fitsread('results/good_seeing/no_dyn_07/contrast.fits');
contrast_dyn = fitsread('results/good_seeing/dyn_07/contrast.fits');



index_LE = 2;
strehl_no_dyn=93.0;
strehl_dyn=93.0;


index_mean = 2;

% avg_con_int_st = mean(contrast_wo_defects(index_mean,4:6));
% avg_con_opt_st = mean(contrast_opt_st(index_mean,4:6));
% avg_con_dd_st = mean(contrast_dd_st(index_mean,4:6));



figure()

% 0.8624
semilogy(contrast_no_dyn(index_mean,:))
hold on;
semilogy(contrast_dyn(index_mean,:))

% semilogy(contrast_three_dead_act_far(index_mean,:))

legend_no_dyn = sprintf('without DM dynamics, strehl = %.1f',strehl_no_dyn);
legend_dyn = sprintf('with DM dynamics, strehl = %.1f',strehl_dyn);

% legend_three_dead_act_far = sprintf('all defects, thre dead act far, strehl = %.1f',strehl_three_dead_act_far);

legend(legend_no_dyn,legend_dyn)
% legend(legend_wo_defects,legend_one_dead_act,legend_two_dead_act_close,legend_two_dead_act_far,legend_three_dead_act_far)

xlabel('r (lamb/D)')
ylabel('contrast')
title('Avg contrast after perfect coronograph')
subtitle("4 kHz, 2 frame delay, wind speed = 6 m/s, seeing = 0.5"", gain = 0.7")
grid()
set(gcf,'Position',[100 100 800 500])
make_it_nicer()