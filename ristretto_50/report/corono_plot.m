corono = fitsread('../results/no_atm/two_dead_act_far/corono.fits');
cor_x = fitsread('../results/no_atm/two_dead_act_far/coroX.fits');
clims = [-7 -4];
figure()
imagesc(cor_x,cor_x,log10(corono),clims)
a = colorbar;
xlabel('x (lamb/D)')
ylabel('y (lamb/D)')
a.Label.String = 'contrast (log scale)';
make_it_nicer()
title('PSF contrast after perfect coronograph')
subtitle("two dead act far")