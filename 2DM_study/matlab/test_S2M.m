% S2M_DM0 = fitsread('../data5/S2M_DM0.fits');
DM0_open = fitsread('../data5/DM0_open.fits');
DM1_open = fitsread('../data5/DM1_open.fits');
DM1_closed = fitsread('../data5/DM1_close.fits');
DM0_closed = fitsread('../data5/DM0_close.fits');
DM1_closed_all = fitsread('../data5/DM1_close_all.fits');
DM0_closed_all = fitsread('../data5/DM0_close_all.fits');

phase = fitsread('../data5/phase.fits');

% S2S = S2M_DM0'*S2M_DM0;
% M2M = S2M_DM0*S2M_DM0';
% M2M = M2M./max(M2M);
% figure()
% imagesc(S2S)
% 
% figure()
% imagesc(M2M)

figure()
plot(abs(DM0_open))
hold on
plot(abs(DM1_open(1:88)))

figure()
plot(abs(DM0_closed))
hold on
plot(abs(DM1_closed(1:88)))

figure()
plot(abs(DM0_closed_all))
hold on
plot(abs(DM1_closed_all(1:88)))

%%
v = 12;
r0 = 0.15;
tau = 0.002;
sigma_time = 6.88*(v*tau/r0).^(5/3);

N = 80;
D = 8;
sigma_fit = 0.335*(D/r0).^(5/3)*N.^(-5/6);

%%
1
0.2             
0.1 0.05

0.2
0.1 0.05

0.2
0.05 0.02

0.2
0.2 0.15
