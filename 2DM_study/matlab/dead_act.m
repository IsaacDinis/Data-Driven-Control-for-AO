act_dead = fitsread('../data4/act_dead.fits');
act_alive = fitsread('../data4/act_alive.fits');
act_dead_1um = fitsread('../data4/act_dead_1um.fits');

%%
fs = 4000;
t = 0:1/fs:10-1/fs;
%%
figure()
plot(t,act_alive)

hold on
plot(t,act_dead)
plot(t,act_dead_1um)

title('Dead actuator neighbour stroke')
legend('actuator alive','dead at 0 um','dead at -1 um','Interpreter','latex','location','northeast');
ylabel('Stroke (um)')
xlabel('Time (s)')
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')