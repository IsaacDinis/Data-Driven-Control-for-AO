contrast_wo_defects = fitsread('results/good_seeing/without_defects/contrast.fits');
contrast_one_dead_act = fitsread('results/good_seeing/one_dead_act/contrast.fits');
contrast_two_dead_act_close = fitsread('results/good_seeing/two_dead_act_close/contrast.fits');
contrast_two_dead_act_far = fitsread('results/good_seeing/two_dead_act_far/contrast.fits');



index_LE = 2;
strehl_wo_defects=93.1;
strehl_one_dead_act=92.1;
strehl_two_dead_act_close=92.0;
strehl_two_dead_act_far=91.9;

index_mean = 2;

% avg_con_int_st = mean(contrast_wo_defects(index_mean,4:6));
% avg_con_opt_st = mean(contrast_opt_st(index_mean,4:6));
% avg_con_dd_st = mean(contrast_dd_st(index_mean,4:6));



figure()

% 0.8624
semilogy(contrast_wo_defects(index_mean,:))
hold on;
semilogy(contrast_one_dead_act(index_mean,:))
semilogy(contrast_two_dead_act_close(index_mean,:))
semilogy(contrast_two_dead_act_far(index_mean,:))


legend_wo_defects = sprintf('without any defects, strehl = %.1f',strehl_wo_defects);
legend_one_dead_act = sprintf('all defects, one dead act, strehl = %.1f',strehl_one_dead_act);
legend_two_dead_act_close = sprintf('all defects, two dead act close, strehl = %.1f',strehl_two_dead_act_close);
legend_two_dead_act_far = sprintf('all defects, two dead act far, strehl = %.1f',strehl_two_dead_act_far);



legend(legend_wo_defects,legend_one_dead_act,legend_two_dead_act_close,legend_two_dead_act_far)

xlabel('r (lamb/D)')
ylabel('contrast')
title('Contrast after perfect coronograph')
subtitle("4 kHz, 1 frame delay, wind speed = 6 m/s, seeing = 0.5""")
grid()
set(gcf,'Position',[100 100 800 500])
make_it_nicer()