
dcao = fitsread('../data/dcao.fits');
dcao_delay_offset = fitsread('../data/dcao_offset.fits');
standalone = fitsread('../data/standalone.fits');
%%
figure()
semilogy(dcao(2,:))
hold on
semilogy(dcao_delay_offset(2,:))
hold on
semilogy(standalone(2,:))

