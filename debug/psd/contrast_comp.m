
dcao = fitsread('../data/dcao.fits');
dcao_delay_offset = fitsread('../data/saxoplus.fits');
standalone = fitsread('../data/standalone.fits');
%%
figure()
semilogy(dcao(2,:))
hold on
semilogy(dcao_delay_offset(2,:))
semilogy(standalone(2,:))

legend('dcao','dcao delay offset','standalone','Interpreter','latex');