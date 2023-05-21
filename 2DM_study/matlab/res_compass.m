res_DM0_alone = fitsread('../data/res_DM0_alone.fits');
res_DM1_alone = fitsread('../data/res_DM1_alone.fits');
res_DM0 = fitsread('../data/res_DM0.fits');
res_DM1 = fitsread('../data/res_DM1.fits');
res_DM0_proj = fitsread('../data/res_DM0_proj.fits');
res_DM0_alone_1kHz = fitsread('../data/res_DM0_alone_1kHz.fits');
dist = fitsread('../data/dist.fits');
dist_DM0 = fitsread('../data/dist_DM0.fits');
dist_DM0_1kHz = fitsread('../data/dist_DM0_1kHz.fits');

t = 0:1/4000:0.25-1/4000;
t_start = 1001;
t_end = 2000;

% stroke_alone = 0.47642;
% stroke = 0.33207;
% strehl_alone = 0.91288 ;
% strehl = 0.93293 ;

%% Dist
figure()
plot(t,dist(t_start:t_end)-mean(dist(t_start:t_end)))
hold on;
ylabel('Dist. amp.')
xlabel('Time (s)')
title('Disturbance')
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/dist.pdf -transparent
%% Intersmapling
figure()

plot(t,res_DM0_alone(t_start:t_end))
hold on
plot(t,res_DM0_alone_1kHz(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Residual')
legend('Res. 1kHz integrator','Res. 4kHz integrator','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/res.pdf -transparent
%% RES alone
figure()

plot(t,res_DM0_alone(t_start:t_end))
hold on
plot(t,res_DM1_alone(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Residual')
legend('Res. 1kHz integrator','Res. 4kHz integrator','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/res.pdf -transparent
%% RES projeted
figure()

plot(t,res_DM0_proj(t_start:t_end))
hold on
plot(t,res_DM1_alone(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Residual')
legend('Res. 1kHz integrator','Res. 4kHz integrator','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/res.pdf -transparent

%% Dist DM0
figure()

plot(t,dist(t_start:t_end)-mean(dist(t_start:t_end)))
hold on
plot(t,dist_DM0(t_start:t_end)-mean(dist_DM0(t_start:t_end)))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Residual')
legend('Res. 1kHz integrator','Res. 4kHz integrator','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/res.pdf -transparent
