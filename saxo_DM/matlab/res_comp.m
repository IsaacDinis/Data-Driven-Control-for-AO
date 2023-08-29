res_cascaded_standalone = load('../data/res_cascaded_standalone.mat').res_cascaded_standalone;
res_cascaded_integrated = load('../data/res_cascaded_integrated.mat').res_cascaded_integrated;
res_1dm = load('../data/res_1dm.mat').res_1dm;
res_1dm_1kHz = load('../data/res_1dm_1kHz.mat').res_1dm;
res_1dm_1kHz_noalias = load('../data/res_1dm_1kHz_noaliasmat.mat').data;

res_parallel_LP = load('../data/res_parallel_LP.mat').data;
res_parallel = load('../data/res_parallel_fixed.mat').data;
res_parallel_noP = load('../data/res_parallel_noP.mat').data;
dist = load('../data/single_mode_dist_nonoise_1_4000_ts.mat').plop.Data;
fs = 4000;
t = 0:1/fs:0.25-1/fs;
t_start = 9001;
t_end = 10000;
%% Dist
figure()
yyaxis left
plot(t,dist(t_start:t_end))
hold on;
ylabel('Dist. amp.')
yyaxis right
plot(t,res_1dm(2,t_start:t_end))
plot(t,res_1dm_1kHz(2,t_start:t_end),color = '#EDB120',LineStyle='-')
xlabel('Time (s)')
ylabel('Res. amp.')
title('Disturbance and residual reference')
legend('Dist.','Res. ref.','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/dist.pdf -transparent
%%
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
%% REs
figure()

plot(t,res_1dm(2,t_start:t_end))
hold on
plot(t,res_1dm_1kHz_noalias(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Residual')
legend('Res. 1kHz integrator','Res. 4kHz integrator','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/res.pdf -transparent

%% Intersampling signal
figure()

plot(t(1:201),res_1dm_1kHz(2,t_start:t_start+200))
hold on
plot(t(1:201),res_1dm_1kHz_noalias(2,t_start+9:t_start+200+9))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Intersampling artifact')
legend('1kHz sampling','4kHz sampling','Interpreter','latex','Location','northwest');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/inter.pdf -transparent

%% Case 1
rms_case1 = rms(res_cascaded_standalone(2,:));
rms_ref = rms(res_1dm(2,:));

figure()
plot(t,res_cascaded_standalone(2,t_start:t_end))
hold on;
plot(t,res_1dm(2,t_start:t_end),'LineWidth',1)
xlabel('Time (s)')
ylabel('Res. amp. (AU)')
title('Residual comparison')
legend('CAO','dCAO','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/case_1.pdf -transparent
%% PSD
n = 40;
w = 1000;

[psd,f] = compute_psd((res_1dm(2,:)-mean(res_1dm(2,:)))',n,w,fs);
sum_psd = sum(psd);
[inter_psd,f] = compute_psd((res_cascaded_standalone(2,:))',n,w,fs);
sum_inter_psd = sum(inter_psd);
figure()
plot(f(1:end),inter_psd(1:end))
hold on;
plot(f(1:end),psd(1:end))
legend('CAO','dCAO','Interpreter','latex');
title('Residual PSD comparison')
xlabel('Frequency (Hz)')
ylabel('Res. amp. (AU)')
% xlim([0,2050])
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
% export_fig ../plot/case_1_psd.pdf -transparent
%% Case 2
figure()
plot(t,res_cascaded_integrated(2,t_start:t_end),'LineWidth',1)
hold on;
plot(t,res_1dm(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Case 2')
legend('Case 2','Ref.','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/case_2.pdf -transparent


%% Case 3
figure()
subplot(1,2,1)
plot(t,res_parallel_noP(2,t_start:t_end))
hold on;
plot(t,res_1dm(2,t_start:t_end))

xlabel('Time (s)')
ylabel('Res. amp.')
legend('Case 3','Ref.','Interpreter','latex');
make_it_nicer()

subplot(1,2,2)
plot(t,res_parallel_noP(2,t_start:t_end))
hold on;

plot(t,res_1dm_1kHz_noalias(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
sgtitle('Case 3','Interpreter','latex','Fontsize',20)
legend('Case 3','1kHz integrator.','Interpreter','latex');
make_it_nicer()

set(gcf, 'Position',  [100, 100, 1400, 450])
set(gcf,'PaperType','A4')

% export_fig ../plot/case_3.pdf -transparent

%% Case 4
figure()
plot(t,res_parallel(2,t_start:t_end),'LineWidth',1)
hold on;
plot(t,res_1dm(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Case 4')
legend('Case 4','Ref.','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/case_4.pdf -transparent

%% Case 5
figure()
plot(t,res_parallel_LP(2,t_start:t_end),'LineWidth',1)
hold on;
plot(t,res_1dm(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Case 5')
legend('Case 5','Ref.','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/case_5.pdf -transparent