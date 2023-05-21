int = load('../data/intersampling.mat').data;
sine = load('../data/intersampling_sine.mat').data;
ramp = load('../data/intersampling_ramp.mat').data;
sine_1234 = load('../data/intersampling_sin_1234.mat').data;
ramp_1234 = load('../data/ramp_1234.mat').data;

ts = 4000;
t = 0:1/ts:0.25-1/ts;
t_start = 9001;
t_end = 10000;

%%
figure()

plot(0:1/ts:0.01-1/ts,int(2,1:40))
% hold on
% plot(t,sine(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
title('Intersampling signal with ramp input')
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
export_fig ../plot/inter_ramp.pdf -transparent

%%
figure()
subplot(1,3,1)
plot(0:1/ts:10,sine(2,:))
% hold on
% plot(t,sine(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')

make_it_nicer()

subplot(1,3,2)
plot(0:1/ts:0.01-1/ts,sine(2,1:40))
% hold on
% plot(t,sine(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
make_it_nicer()

subplot(1,3,3)
plot(1.56:1/ts:1.56+0.015-1/ts,sine(2,6251:6310))
% hold on
% plot(t,sine(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
make_it_nicer()

sgtitle('Intersampling signal with sine input','Interpreter','latex','Fontsize',20)

set(gcf, 'Position',  [100, 100, 2100, 450])
set(gcf,'PaperType','A4')
export_fig ../plot/inter_sine.pdf -transparent -nocrop

%%
figure()

plot(0:1/ts:0.01-1/ts,ramp_1234(2,1:40))
% hold on
% plot(t,ramp_1234(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
title('Intersampling signal with ramp input')
make_it_nicer()

set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')

export_fig ../plot/inter_ramp_1234.pdf -transparent