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




% legend_one_dead_act = sprintf('one dead act, strehl = %.2f',strehl_int_st);
% legend_two_dead_act_close = sprintf('two dead act close, strehl = %.2f',strehl_opt_st);
% legend_two_dead_act_far = sprintf('two dead act far, strehl = %.2f',strehl_dd_st);
% legend_without_defects = sprintf('without deffects, strehl = %.2f',strehl_int_dcao);

% legend_one_dead_act = sprintf('one dead act');
% legend_two_dead_act_close = sprintf('two dead act close');
% legend_two_dead_act_far = sprintf('two dead act far');
% legend_without_defects = sprintf('without defects');

legend_one_dead_act = sprintf('one dead act');
legend_two_dead_act_close = sprintf('two dead act close');
legend_two_dead_act_far = sprintf('two dead act far');
legend_without_defects = sprintf('without defects');



legend(legend_without_defects,legend_one_dead_act,legend_two_dead_act_close,legend_two_dead_act_far)

xlabel('r (lamb/D)')
ylabel('contrast')
title('Radial averaged contrast after perfect coronograph')
subtitle("lambda = 750 nm, wind speed = 6 m/s, seeing = 0.5""")
grid()
set(gcf,'Position',[100 100 900 600])
make_it_nicer()