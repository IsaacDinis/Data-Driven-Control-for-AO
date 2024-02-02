x = [-3:.1:3];
lodm = normpdf(x,0,1);
hodm = normpdf(x,0,1/2.2);
hodm_l = normpdf(x,-1,1/2.2);
hodm_r = normpdf(x,1,1/2.2);
hodm = hodm/max(hodm)*max(lodm)*0.9;
hodm_l = hodm_l/max(hodm_l)*lodm(x==-1)*0.9;
hodm_r = hodm_r/max(hodm_r)*lodm(x==1)*0.9;
%%
figure()
plot(x,lodm)
hold on;
plot(x,hodm)
plot(x,hodm_l)
plot(x,hodm_r)

figure()
plot(x,lodm)
hold on;
plot(x,hodm+hodm_l+hodm_r)

figure()
plot(-lodm+hodm_l+hodm_r)

figure()
plot(-lodm+hodm_l+hodm_r+hodm)

%%
figure()
subplot(1,3,1)
plot(x,hodm)
ylim([-0.5,0.5])
legend('Goal to achieve')

subplot(1,3,2)
plot(x,lodm)
hold on;
plot(x,-hodm_l)
plot(x,-hodm_r)
ylim([-0.5,0.5])
legend('LODM','left neighbour','right neighbour')

subplot(1,3,3)
plot(lodm-hodm_l-hodm_r)
ylim([-0.5,0.5])
legend(['Result'])

set(gcf,'Position',[100 100 750 250])
sgtitle('Phantom actuator')
%%
figure()
subplot(1,4,1)
plot(x,hodm)
ylim([-0.5,0.5])
legend('Dead actuator')

subplot(1,4,2)
plot(x,-lodm)
hold on;
plot(x,hodm_l)
plot(x,hodm_r)
ylim([-0.5,0.5])
legend('LODM','left neighbour','right neighbour')

subplot(1,4,3)
plot(-lodm+hodm_l+hodm_r)
ylim([-0.5,0.5])
legend(['Sum of LODM, left' newline 'and right neighbour'])

subplot(1,4,4)
plot((-lodm+hodm_l+hodm_r)+hodm)
legend('Sum of everything')
ylim([-0.5,0.5])
set(gcf,'Position',[100 100 1000 250])
sgtitle('Achieve flat with a dead actuator')
%%

figure()
subplot(1,2,1)
plot(x,-lodm*0.5)
hold on;
plot(x,hodm)
plot(x,hodm_l*0.5)
plot(x,hodm_r*0.5)
ylim([-0.5,0.5])
legend('LODM','dead actuator', 'left neighbour','right neighbour')

subplot(1,2,2)
plot((-lodm+hodm_l+hodm_r)*0.5+hodm)
ylim([-0.5,0.5])
legend('Sum of everything')

set(gcf,'Position',[100 100 500 250])
sgtitle('Decrease dead actuator stroke')
%%

figure()
subplot(1,2,1)
plot(x,-lodm*1.5)
hold on;
plot(x,hodm)
plot(x,hodm_l*1.5)
plot(x,hodm_r*1.5)
ylim([-0.8,0.8])
legend('LODM','dead actuator', 'left neighbour','right neighbour')

subplot(1,2,2)
plot((-lodm+hodm_l+hodm_r)*1.5+hodm)
legend('Sum of everything')
ylim([-0.8,0.8])
set(gcf,'Position',[100 100 500 250])
sgtitle('Inverse dead actuator stroke')
%%

figure()
subplot(1,2,1)
plot(x,-lodm*1.5)
hold on;
plot(x,hodm)
plot(x,hodm_l*2)
plot(x,hodm_r*2)
ylim([-0.8,0.8])
legend('LODM','dead actuator', 'left neighbour','right neighbour')
subplot(1,2,2)
plot((-lodm+hodm_l+hodm_r)*2+hodm)
ylim([-0.8,0.8])
set(gcf,'Position',[100 100 500 250])
legend('Sum of everything')
sgtitle('Impossible combination')