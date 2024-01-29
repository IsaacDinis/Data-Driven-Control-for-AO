saxoplus = fitsread('../results/bright1_1_3ms/dcao/contrastCurves/saxoplus.fits');
saxo = fitsread('../results/bright1_1_3ms/dcao/contrastCurves/saxo.fits');
saxoplus2 = fitsread('../results/bright1_1_3ms/dcao_dd_full/contrastCurves/saxoplus.fits');

figure()
semilogy(saxo(2,:))
hold on;
semilogy(saxoplus2(2,:))
semilogy(saxoplus(2,:))
legend('saxo alone','data-driven','integrator')
xlabel('r (lamb/D)')
ylabel('contrast')
title('contrast after perfect coronograph')