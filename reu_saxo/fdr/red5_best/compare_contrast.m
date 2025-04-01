contrast_ref_no_vib = fitsread('results/ref_no_vib/contrastCurves/saxoplus.fits');
contrast_ref_vib = fitsread('results/ref_vib/contrastCurves/saxoplus.fits');
contrast_dd_no_vib = fitsread('results/dd4ao_no_vib/contrastCurves/saxoplus.fits');
contrast_dd_vib = fitsread('results/dd4ao_vib/contrastCurves/saxoplus.fits');

strehl_ref_no_vib = fitsread('results/ref_no_vib/strehlCurves/saxoplus.fits');
strehl_ref_vib = fitsread('results/ref_vib/strehlCurves/saxoplus.fits');
strehl_dd_no_vib = fitsread('results/dd4ao_no_vib/strehlCurves/saxoplus.fits');
strehl_dd_vib = fitsread('results/dd4ao_vib/strehlCurves/saxoplus.fits');

index_LE = 2;
strehl_ref_no_vib=strehl_ref_no_vib(index_LE,end)*100;
strehl_ref_vib=strehl_ref_vib(index_LE,end)*100;
strehl_dd_no_vib=strehl_dd_no_vib(index_LE,end)*100;
strehl_dd_vib=strehl_dd_vib(index_LE,end)*100;

index_mean = 2;

avg_con_ref_no_vib = mean(contrast_ref_no_vib(index_mean,4:6));
avg_con_ref_vib = mean(contrast_ref_vib(index_mean,4:6));
avg_con_dd_no_vib = mean(contrast_dd_no_vib(index_mean,4:6));
avg_con_dd_vib = mean(contrast_dd_vib(index_mean,4:6));


figure()
index_mean = 2;
% 0.8624
end_contrast = 20;
semilogy(contrast_ref_no_vib(index_mean,1:end_contrast),'r')
hold on;
semilogy(contrast_ref_vib(index_mean,1:end_contrast),'r--')
semilogy(contrast_dd_no_vib(index_mean,1:end_contrast),'b')
semilogy(contrast_dd_vib(index_mean,1:end_contrast),'b--')


legend_ref_no_vib = sprintf('ref controller, strehl = %.2f',strehl_ref_no_vib);
legend_ref_vib = sprintf('with vibrations, strehl = %.2f',strehl_ref_vib);
legend_dd_no_vib = sprintf('dd4ao controller, strehl = %.2f',strehl_dd_no_vib);
legend_dd_vib = sprintf('with vibrations, strehl = %.2f',strehl_dd_vib);



legend(legend_ref_no_vib,legend_ref_vib,legend_dd_no_vib,legend_dd_vib)

xlabel('r (lamb/D)')
ylabel('contrast')
title('Red 5 best, saxo+ contrast after perfect coronograph')
subtitle("G mag = 16.8, J mag = 12.5, t0 = 9 ms, seeing = 0.4""")
grid()
set(gcf,'Position',[100 100 1100 700])
make_it_nicer()