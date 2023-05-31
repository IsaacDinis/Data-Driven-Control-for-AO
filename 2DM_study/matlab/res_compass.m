res_DM0_alone = fitsread('../data/res_DM0_alone.fits');
res_DM1_alone = fitsread('../data/res_DM1_alone.fits');
res_DM0 = fitsread('../data/res_DM0.fits');
res_DM1 = fitsread('../data/res_DM1.fits');
res_DM0_proj = fitsread('../data/res_DM0_proj.fits');
res_DM0_alone_1kHz = fitsread('../data/res_DM0_alone_1kHz.fits');
dist = fitsread('../data/dist.fits');
dist_DM0 = fitsread('../data/dist_DM0.fits');
dist_DM0_1kHz = fitsread('../data/dist_DM0_1kHz.fits');
% res_DM1_alone_proj_1kHz = fitsread('../data/res_DM1_alone_projeted_1kHz.fits');
res_DM1_alone_proj = fitsread('../data/res_DM1_alone_projeted.fits');
t = 0:1/4000:0.25-1/4000;
t_start = 1001;
t_end = 2000;

% stroke_alone = 0.47642;
% stroke = 0.33207;
% strehl_alone = 0.91288 ;
% strehl = 0.93293 ;
% strehl LODM = 0.62395 
%% Dist
figure()
plot(0:1/4000:8-1/4000,dist)
hold on;
plot(0:1/4000:8-1/4000,dist_DM0)
ylabel('Dist. amp.')
xlabel('Time (s)')
title('Disturbance projected on DMs subspaces')
legend('HODM 1st KL','LODM 1st KL','Interpreter','latex');
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
figure()
plot(t,res_DM0_alone_1kHz(t_start:t_end),'LineWidth',1)
hold on
plot(t,res_DM0_alone(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Residual on LODM 1st KL')
legend('1 kHz sampling','4 kHz sampling','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/res.pdf -transparent
%% PSD
n = 20;
w = 100;

[psd,f] = compute_psd((res_DM0_alone_1kHz)',n,w,fs);

[inter_psd,f] = compute_psd((res_DM0_alone)',n,w,fs);

figure()
plot(f(1:end),psd(1:end))
hold on;
plot(f(1:end),inter_psd(1:end))

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
figure()
subplot(1,2,1)
plot(t,res_DM1(t_start:t_end))
hold on
plot(t,res_DM1_alone(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
legend('HODM only','Both DMs','Interpreter','latex');
make_it_nicer()
subplot(1,2,2)
plot(t,res_DM0_proj(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
legend('LODM only','Interpreter','latex');
sgtitle('Residual on HODM 1st KL','Interpreter','latex','Fontsize',20)
set(gcf, 'Position',  [100, 100, 1400, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/res_DM1.pdf -transparent -nocrop
%% RES tot DM0
figure()

plot(t,res_DM0_alone_1kHz(t_start:t_end))
hold on
plot(t,res_DM1_alone_proj(t_start:t_end))
plot(t,res_DM0(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Residual on LODM 1st KL')
legend('LODM only','HODM only','both DMs','Interpreter','latex', 'Location','southwest');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/compass/res_DM0.pdf -transparent