dist_DM1 = fitsread('../data_parallel/dist_DM1.fits');
dist_DM0 = fitsread('../data_parallel/dist_DM0.fits');
dist_tilt = fitsread('../data_parallel/dist_tilt.fits');

factor = max(abs(dist_tilt))/max(abs(dist_DM0));
dist_DM0 = dist_DM0*factor;
dist_DM1 = dist_DM1*factor;

res_DM0_alone = fitsread('../data_parallel/res_DM0_alone.fits')*factor;
res_DM1_alone = fitsread('../data_parallel/res_DM1_alone.fits')*factor;

res_DM0 = fitsread('../data_parallel/res_DM0.fits')*factor;
res_DM1 = fitsread('../data_parallel/res_DM1.fits')*factor;

res_tilt = fitsread('../data_parallel/res_tilt.fits');
res_tilt_DM1 = fitsread('../data_parallel/res_tilt_DM1.fits');
res_tilt_DM0 = fitsread('../data_parallel/res_tilt_DM0.fits');
% res_tilt_DM0_4kHz = -fitsread('../data_parallel/res_tilt_DM0.fits');

res_DM0_proj = fitsread('../data_parallel/res_DM0_proj.fits')*factor;
res_DM1_proj = fitsread('../data_parallel/res_DM1_proj.fits')*factor;

command_DM0 = fitsread('../data_parallel/command_DM0.fits')*factor;
command_DM1 = fitsread('../data_parallel/command_DM1.fits')*factor;
command_DM0_alone = fitsread('../data_parallel/command_DM0_alone.fits')*factor;
command_DM1_alone = fitsread('../data_parallel/command_DM1_alone.fits')*factor;
command_DM1_alone = command_DM1_alone(1001:end);
command_DM1 = command_DM1(1001:end);
% 0.02309
% 0.00006 
% res_DM1_integrating_DM0 = fitsread('../data2/res_DM1_integrating_DM0.fits');
% M2M = fitsread('../data2/M2M.fits');


% dist_DM0_1kHz = fitsread('../data2/dist_DM0_1kHz.fits');
% res_DM1_alone_proj_1kHz = fitsread('../data/res_DM1_alone_projeted_1kHz.fits');

t = 0:1/4000:1-1/4000;
t_start = 1001;
t_end = 2000;
start = 1;
% stroke_alone = 0.47874 ;
% stroke = 0.32459 
% strehl_alone = 0.86682;
% strehl =  0.93649;
% strehl LODM = 0.62814 
%% Dist
figure()
plot(t(start:end),dist_DM1(start:end))
hold on;
plot(t(start:end),dist_DM0(start:end))
plot(t(start:end),dist_tilt(start:end))
ylabel('Dist. amp. (nm)')
xlabel('Time (s)')
title('Disturbance')
legend('projeted on slow TT','projeted on fast TT','phase tilt','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/dist.pdf -transparent





%% RES tot DM1

figure()

plot(t(start:end),res_DM0_proj(start:end))
hold on
plot(t(start:end),res_DM1_alone(start:end))
plot(t(start:end),res_DM1(start:end))


xlabel('Time (s)')
ylabel('Res. amp. (nm)')
legend('slow TT only','fast TT only','both TTs','Interpreter','latex');
make_it_nicer()


title('Residual on fast TT mirror','Interpreter','latex','Fontsize',20)
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/res_DM1.pdf -transparent -nocrop
%% RES tot DM0
% figure()
% 
% plot(t(20:end),res_DM0_alone(20:end))
% hold on
% plot(t(20:end),res_DM1_proj(20:end))
% plot(t(20:end),res_DM0(20:end))
% % plot(t,res_DM1_integrating_DM0(t_start:t_end))
% xlabel('Time (s)')
% ylabel('Res. amp. (nm)')
% title('Residual on LODM 1st KL')
% legend('LODM only','HODM only','both DMs','Interpreter','latex');
% make_it_nicer()
% set(gcf, 'Position',  [100, 100, 700, 450])
% set(gcf,'PaperType','A4')
% make_it_nicer()
% % export_fig ../plot/compass/res_DM0.pdf -transparent

%% RES tot DM0

figure()

plot(t(start:end),res_DM0_alone(start:end))
hold on
plot(t(start:end),res_DM1_proj(start:end))
plot(t(start:end),res_DM0(start:end))


xlabel('Time (s)')
ylabel('Res. amp. (nm)')
legend('slow TT only','fast TT only','both TTs','Interpreter','latex');
make_it_nicer()


title('Residual on slow TT mirror','Interpreter','latex','Fontsize',20)
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/res_DM0.pdf -transparent -nocrop
%% M2M
% figure()
% plot(M2M(:,1))
% title('LODM 1 KL with HODM KLs')
% xlabel('HODM KL mode')
% make_it_nicer()
% set(gcf, 'Position',  [100, 100, 700, 450])
% set(gcf,'PaperType','A4')
% % export_fig ../plot/compass/M2M.pdf -transparent




%% RES tot tilt
figure()
plot(t(start:end),res_tilt_DM0(start:end))
hold on
plot(t(start:end),res_tilt_DM1(start:end))
plot(t(start:end),res_tilt(start:end))
% plot(t,res_DM1_integrating_DM0(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp. (um)')
title('Residual on tilt')
legend('slow TT only','fast TT only','both TTs','Interpreter','latex');
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

legend('slow TT only','fast TT only','both TTs','Interpreter','latex');
title('Tilt residual PSD')
xlabel('Frequency (Hz)')
ylabel('Amp.')
% xlim([0,2050])
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
% export_fig ../plot/compass/tilt_psd.pdf -transparent
%% COMMAND
figure()

plot(t(start:end),command_DM0(start:end))
hold on
% plot(t(start:end),command_DM0_alone(start:end))
plot(t(start:end),command_DM1(start:end))
plot(t(start:end),command_DM1_alone(start:end))


xlabel('Time (s)')
ylabel('Res. amp. (nm)')
% legend('DM0','DM1','DM1 alone','Interpreter','latex');
make_it_nicer()


title('Residual on fast TT mirror','Interpreter','latex','Fontsize',20)
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()