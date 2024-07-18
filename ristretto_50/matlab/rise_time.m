ol_dyn = fitsread('ol_dyn.fits');
ol_no_dyn = fitsread('ol_no_dyn.fits');
cl_dyn = fitsread('cl_dyn.fits');
cl_no_dyn = fitsread('cl_no_dyn.fits');
ol_dyn = [ol_dyn,ones(1,15)*ol_dyn(end)];
ol_no_dyn = [ol_no_dyn,ones(1,15)*ol_no_dyn(end)];
fs = 4000;
t = linspace(0,19/fs,20);

figure()
plot(t,ol_dyn);

hold on;

plot(t,cl_dyn);
plot(t,cl_no_dyn);
plot(t,ol_no_dyn);
legend('open loop with DM dynamics','closed loop with DM dynamics','closed loop without DM dynamics','open loop without DM dynamics')
xlabel('time (s)')
title('tilt step response, compass')
make_it_nicer()