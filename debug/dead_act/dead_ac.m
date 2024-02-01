x = [-3:.1:3];
lodm = normpdf(x,0,1);
hodm = normpdf(x,0,1/2.2);
hodm_l = normpdf(x,-1,1/2.2);
hodm_r = normpdf(x,1,1/2.2);
hodm = hodm/max(hodm)*max(lodm)*0.9;
hodm_l = hodm_l/max(hodm_l)*lodm(x==-1)*0.9;
hodm_r = hodm_r/max(hodm_r)*lodm(x==1)*0.9;
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

figure()
subplot(1,4,1)
plot(x,hodm)
ylim([-0.5,0.5])

subplot(1,4,2)
plot(x,-lodm)
hold on;
plot(x,hodm_l)
plot(x,hodm_r)
ylim([-0.5,0.5])

subplot(1,4,3)
plot(-lodm+hodm_l+hodm_r)
ylim([-0.5,0.5])

subplot(1,4,4)
plot((-lodm+hodm_l+hodm_r)-0.5+hodm)
ylim([-0.5,0.5])

figure()
subplot(1,2,1)
plot(x,-lodm*0.5)
hold on;
plot(x,hodm)
plot(x,hodm_l*0.5)
plot(x,hodm_r*0.5)
ylim([-0.5,0.5])

subplot(1,2,2)
plot((-lodm+hodm_l+hodm_r)*0.5+hodm)
ylim([-0.5,0.5])

figure()
subplot(1,2,1)
plot(x,-lodm*1.5)
hold on;
plot(x,hodm)
plot(x,hodm_l*1.5)
plot(x,hodm_r*1.5)
ylim([-0.5,0.5])

subplot(1,2,2)
plot((-lodm+hodm_l+hodm_r)*1.5+hodm)
ylim([-0.5,0.5])

figure()
subplot(1,2,1)
plot(x,-lodm*1.5)
hold on;
plot(x,hodm)
plot(x,hodm_l*2)
plot(x,hodm_r*2)
ylim([-0.5,0.5])

subplot(1,2,2)
plot((-lodm+hodm_l+hodm_r)*2+hodm)
ylim([-0.5,0.5])
