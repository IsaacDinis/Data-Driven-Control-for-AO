res_cascaded_standalone = load('../data/res_cascaded_standalone.mat').res_cascaded_standalone;
res_cascaded_integrated = load('../data/res_cascaded_integrated.mat').res_cascaded_integrated;
res_1dm = load('../data/res_1dm.mat').res_1dm;
res_1dm_1kHz = load('../data/res_1dm_1kHz.mat').res_1dm;

% res_parallel_LP = load('../data/res_parallel_LP.mat').data;
% res_parallel = load('../data/res_parallel.mat').data;
% res_parallel_noP = load('../data/res_parallel_noP.mat').data;
dist = load('../data/single_mode_dist_nonoise_1_4000_ts.mat').plop.Data;

t = 0:1/4000:0.25-1/4000;
%% Dist
figure()
yyaxis left
plot(t,dist(9001:10000))
hold on;
ylabel('Dist. amp.')
yyaxis right
plot(t,res_1dm(2,9001:10000))
plot(t,res_1dm_1kHz(2,9001:10000),color = '#EDB120',LineStyle='-')
xlabel('Time (s)')
ylabel('Res. amp.')
title('Disturbance and residual reference')
legend('Dist.','Res. ref.','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/dist.pdf -transparent

%% Case 1
figure()
plot(t,res_cascaded_standalone(2,9001:10000))
hold on;
plot(t,res_1dm(2,9001:10000),'LineWidth',1)
xlabel('Time (s)')
ylabel('Res. amp.')
title('Case 1')
legend('Case 1','Ref.','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/case_1.pdf -transparent

%% Case 2
figure()
plot(t,res_cascaded_integrated(2,9001:10000),'LineWidth',1)
hold on;
plot(t,res_1dm(2,9001:10000))
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
plot(t,res_parallel_noP(2,9001:10000))
hold on;
plot(t,res_1dm(2,9001:10000))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Case 3')
legend('Case 3','Ref.','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/case_3.pdf -transparent

%% Case 4
figure()
plot(t,res_parallel(2,9001:10000),'LineWidth',1)
hold on;
plot(t,res_1dm(2,9001:10000))
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
plot(t,res_parallel_LP(2,9001:10000),'LineWidth',1)
hold on;
plot(t,res_1dm(2,9001:10000))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Case 5')
legend('Case 5','Ref.','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/case_5.pdf -transparent