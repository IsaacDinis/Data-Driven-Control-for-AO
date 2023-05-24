int_1 = load('../data/cl_inter_1st.mat').data;
int_2 = load('../data/cl_inter_2st.mat').data;

int_1 = int_1(2,:);
int_2 = int_2(2,:);

fs = 4000;
t = 0:1/fs:0.25-1/fs;
t_start = 9001;
t_end = 10000;

%%
figure()

plot(0:1/fs:0.01-1/fs,int_1(1:40))
% hold on
% plot(t,sine(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
title('intersampling signal with ramp input')
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/int_1er_ramp.pdf -transparent
%%
figure()

plot(0:1/fs:0.01-1/fs,int_2(1:40))
% hold on
% plot(t,sine(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
title('intersampling signal with ramp input')
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/int_1er_ramp.pdf -transparent

%% PSD
n = 40;
w = 1000;

[int_1er_psd_sine,f] = compute_psd((int_1-mean(int_1))',n,w,fs);


figure()
plot(f(1:end),int_1er_psd_sine(1:end))
title('ramp (1000 Hz) input')
xlabel('Frequency (Hz)')
ylabel('Amp.')
xlim([0,2050])
make_it_nicer()
set(gcf, 'Position',  [100, 100, 1500, 800])
set(gcf,'PaperType','A4')

%% PSD
n = 40;
w = 1000;

[int_1er_psd_sine,f] = compute_psd((int_2-mean(int_2))',n,w,fs);


figure()
plot(f(1:end),int_1er_psd_sine(1:end))
title('ramp (1000 Hz) input')
xlabel('Frequency (Hz)')
ylabel('Amp.')
xlim([0,2050])
make_it_nicer()
set(gcf, 'Position',  [100, 100, 1500, 800])
set(gcf,'PaperType','A4')