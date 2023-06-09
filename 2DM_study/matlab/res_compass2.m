res_DM0_alone = fitsread('../data2/res_DM0_alone.fits')*1000;
res_DM1_alone = fitsread('../data2/res_DM1_alone.fits')*1000;

res_DM0 = fitsread('../data2/res_DM0.fits')*1000;
res_DM1 = fitsread('../data2/res_DM1.fits')*1000;

res_tilt = fitsread('../data2/res_tilt.fits')*1000;
res_tilt_DM1 = fitsread('../data2/res_tilt_DM1.fits')*1000;
res_tilt_DM0 = fitsread('../data2/res_tilt_DM0.fits')*1000;
res_tilt_DM0_4kHz = fitsread('../data2/res_tilt_DM0_4kHz.fits')*1000;

res_DM0_proj = fitsread('../data2/res_DM0_proj.fits')*1000;
res_DM1_proj = fitsread('../data2/res_DM1_proj.fits')*1000;

% res_DM1_integrating_DM0 = fitsread('../data2/res_DM1_integrating_DM0.fits');
M2M = fitsread('../data2/M2M.fits');

dist_DM1 = fitsread('../data2/dist_DM1.fits')*1000;
dist_DM0 = fitsread('../data2/dist_DM0.fits')*1000;
dist_tilt = fitsread('../data2/dist_tilt.fits')*1000;

% dist_DM0_1kHz = fitsread('../data2/dist_DM0_1kHz.fits');
% res_DM1_alone_proj_1kHz = fitsread('../data/res_DM1_alone_projeted_1kHz.fits');

t = 0:1/4000:0.25-1/4000;
t_start = 1001;
t_end = 2000;

% stroke_alone = 0.47874 ;
% stroke = 0.32459 
% strehl_alone = 0.86682;
% strehl =  0.93649;
% strehl LODM = 0.62814 
%% Dist
figure()
plot(0:1/4000:1-1/4000,dist_DM1)
hold on;
plot(0:1/4000:1-1/4000,dist_DM0)
plot(0:1/4000:1-1/4000,dist_tilt)
ylabel('Dist. amp. (nm)')
xlabel('Time (s)')
title('Disturbance')
legend('HODM 1st KL','LODM 1st KL','phase tilt','Interpreter','latex','location','northwest');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/dist.pdf -transparent
%% Dist 2
% figure()
% 
% plot(t,dist(t_start:t_end)-mean(dist(t_start:t_end)))
% hold on
% plot(t,dist_DM0(t_start:t_end)-mean(dist_DM0(t_start:t_end)))
% xlabel('Time (s)')
% ylabel('Dist. amp.')
% title('Disturbance')
% legend('Dist. on HODM 1st KL','Dist. on LODM 1st KL','Interpreter','latex', Location='northwest');
% make_it_nicer()
% set(gcf, 'Position',  [100, 100, 700, 450])
% set(gcf,'PaperType','A4')
% make_it_nicer()
% % export_fig ../plot/res.pdf -transparent

%% Intersmapling
% figure()
% plot(t,res_DM0_alone_1kHz(t_start:t_end),'LineWidth',1)
% hold on
% plot(t,res_DM0_alone(t_start:t_end))
% xlabel('Time (s)')
% ylabel('Res. amp.')
% title('Residual on LODM 1st KL')
% legend('1 kHz sampling','4 kHz sampling','Interpreter','latex');
% make_it_nicer()
% set(gcf, 'Position',  [100, 100, 700, 450])
% set(gcf,'PaperType','A4')
% make_it_nicer()
% export_fig ../plot/compass/res.pdf -transparent
% %% PSD
% n = 20;
% w = 100;
% 
% [psd,f] = compute_psd((res_DM0_alone_1kHz)',n,w,fs);
% 
% [inter_psd,f] = compute_psd((res_DM0_alone)',n,w,fs);
% 
% figure()
% plot(f(1:end),psd(1:end))
% hold on;
% plot(f(1:end),inter_psd(1:end))

%% RES alone
% figure()
% 
% plot(t,res_DM0_proj(t_start:t_end))
% xlabel('Time (s)')
% ylabel('Res. amp.')
% title('Residual')
% legend('Res. 1kHz integrator','Res. 4kHz integrator','Interpreter','latex');
% make_it_nicer()
% set(gcf, 'Position',  [100, 100, 700, 450])
% set(gcf,'PaperType','A4')
% make_it_nicer()
% export_fig ../plot/res.pdf -transparent
%% RES projeted
% figure()
% 
% plot(t,res_DM0_proj(t_start:t_end))
% hold on
% plot(t,res_DM1_alone(t_start:t_end))
% xlabel('Time (s)')
% ylabel('Res. amp.')
% title('Residual')
% legend('Res. 1kHz integrator','Res. 4kHz integrator','Interpreter','latex');
% make_it_nicer()
% set(gcf, 'Position',  [100, 100, 700, 450])
% set(gcf,'PaperType','A4')
% make_it_nicer()
% % export_fig ../plot/res.pdf -transparent
%% RES tot DM1
start = 30;
figure()
subplot(1,2,1)
plot(t(start:end),res_DM1(start:end))
hold on
plot(t(start:end),res_DM1_alone(start:end))
xlabel('Time (s)')
ylabel('Res. amp. (nm)')
legend('HODM only','Both DMs','Interpreter','latex');
make_it_nicer()
subplot(1,2,2)
plot(t(start:end),res_DM0_proj(start:end))
xlabel('Time (s)')
ylabel('Res. amp. (nm)')
legend('LODM only','Interpreter','latex');
sgtitle('Residual on HODM 1st KL','Interpreter','latex','Fontsize',20)
set(gcf, 'Position',  [100, 100, 1400, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/res_DM1.pdf -transparent -nocrop
%% RES tot DM0
figure()

plot(t(20:end),res_DM0_alone(20:end))
hold on
plot(t(20:end),res_DM1_proj(20:end))
plot(t(20:end),res_DM0(20:end))
% plot(t,res_DM1_integrating_DM0(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp. (nm)')
title('Residual on LODM 1st KL')
legend('LODM only','HODM only','both DMs','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/res_DM0.pdf -transparent

%% M2M
figure()
plot(M2M(:,1))
title('LODM 1 KL with HODM KLs')
xlabel('HODM KL mode')
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
% export_fig ../plot/compass/M2M.pdf -transparent




%% RES tot DM0
figure()

plot(t,res_tilt_DM0)
hold on
plot(t,res_tilt_DM1)
plot(t,res_tilt)
plot(t,res_tilt_DM0_4kHz)
% plot(t,res_DM1_integrating_DM0(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp. (um)')
title('Residual on tilt')
legend('LODM only','HODM only','both DMs','LODM only 4 kHz','Interpreter','latex', 'Location','southwest');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/res_DM0.pdf -transparent

%% RES tot DM0
t = 0:1/4000:1-1/4000;
figure()

plot(t,res_tilt_DM0)
hold on
plot(t,res_tilt_DM1)
plot(t,res_tilt)
plot(t,res_tilt_DM0_4kHz)
% plot(t,res_DM1_integrating_DM0(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp. (nm)')
title('Tilt residual')
legend('LODM only','HODM only','both DMs','LODM only 4 kHz','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/res_tilt.pdf -transparent
%%
n = 10;
w = 399;
fs = 4000;
[psd_tilt_DM0,f] = compute_psd(res_tilt_DM0',n,w,fs);
[psd_tilt_DM1,f] = compute_psd(res_tilt_DM1',n,w,fs);
[psd_tilt,f] = compute_psd(res_tilt',n,w,fs);


figure()
loglog(f(1:end),psd_tilt_DM0(1:end))
hold on;
loglog(f(1:end),psd_tilt_DM1(1:end))
loglog(f(1:end),psd_tilt(1:end))

legend('LODM only','HODM only','both DMs','Interpreter','latex');
title('Tilt residual PSD')
xlabel('Frequency (Hz)')
ylabel('Amp.')
% xlim([0,2050])
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
export_fig ../plot/tilt_psd.pdf -transparent