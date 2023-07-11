act_dead = fitsread('../data4/act_dead.fits');
act_alive = fitsread('../data4/act_alive.fits');
act_dead_1um = fitsread('../data4/act_dead_1um.fits');
figure()
plot(act_dead)
hold on
plot(act_alive)
plot(act_dead_1um)