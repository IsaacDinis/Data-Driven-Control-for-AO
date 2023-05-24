stroke_cascaded_standalone = load('../data/stroke_cascaded_standalone.mat').data;
stroke_cascaded_integrated = load('../data/stroke_cascaded_integrated.mat').data;
stroke_parallel_standalone = load('../data/stroke_parallel_standalone.mat').data;
stroke_parallel_integrated = load('../data/stroke_parallel_integrated.mat').data;



fs = 4000;
t = 0:1/fs:2.5-1/fs;
t_start = 100;
t_end = 10099;
% %% Dist
% figure()
% yyaxis left
% plot(t,dist(t_start:t_end))
% hold on;
% ylabel('Dist. amp.')
% yyaxis right
% plot(t,stroke_1dm(2,t_start:t_end))
% plot(t,res_1dm_1kHz(2,t_start:t_end),color = '#EDB120',LineStyle='-')
% xlabel('Time (s)')
% ylabel('Res. amp.')
% title('Disturbance and residual reference')
% legend('Dist.','Res. ref.','Interpreter','latex');
% make_it_nicer()
% set(gcf, 'Position',  [100, 100, 700, 450])
% set(gcf,'PaperType','A4')
% make_it_nicer()
% % export_fig ../plot/dist.pdf -transparent
%%
% figure()
% plot(t,dist(t_start:t_end)-mean(dist(t_start:t_end)))
% hold on;
% ylabel('Dist. amp.')
% xlabel('Time (s)')
% title('Disturbance')
% make_it_nicer()
% set(gcf, 'Position',  [100, 100, 700, 450])
% set(gcf,'PaperType','A4')
% make_it_nicer()
% % export_fig ../plot/dist.pdf -transparent
%% Dist
rms_case1 = rms(stroke_cascaded_standalone(2,:));
rms_ref = rms(stroke_1dm(2,:));

figure()
yyaxis left
plot(t,stroke_1dm(2,t_start:t_end))
hold on;
ylabel('Amp.')
yyaxis right
plot(t,stroke_cascaded_standalone(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
title('Stroke')
legend('Ref.','Case 1','Interpreter','latex',Location='northwest');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/stroke_case1.pdf -transparent
%% REs


rms_case2 = rms(stroke_cascaded_integrated(2,:));
rms_ref = rms(stroke_1dm(2,:));

figure()
subplot(1,2,1)
yyaxis left
plot(t,stroke_cascaded_integrated(2,t_start:t_end))
hold on;
ylabel('Amp.')
yyaxis right
plot(t,stroke_1dm(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
legend('Case 2','Ref','Interpreter','latex',Location='northwest');
make_it_nicer()

subplot(1,2,2)
t1 = 0:1/fs:0.05-1/fs;
t_start1 = 9001;
t_end1 = 9200;

plot(t1,stroke_cascaded_integrated(2,t_start1:t_end1))
hold on;
plot(t1,stroke_cascaded_standalone(2,t_start1:t_end1))
xlabel('Time (s)')
ylabel('Amp.')

legend('Case 2','Case 1','Interpreter','latex',Location='northwest');
make_it_nicer()
sgtitle('Stroke','Interpreter','latex','Fontsize',20)
set(gcf, 'Position',  [100, 100, 1400, 450])
set(gcf,'PaperType','A4')

% export_fig ../plot/stroke_case2.pdf -transparent -nocrop

%% Case 1
figure()

plot(t,stroke_cascaded_integrated(2,t_start:t_end))
hold on
plot(t,stroke_cascaded_standalone(2,t_start:t_end))
xlabel('Time (s)')
ylabel('Amp.')
title('Stroke')
legend('Ref.','Case 1','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/res.pdf -transparent

%% Case 2
rms_case1 = rms(stroke_cascaded_standalone(2,:));
rms_ref = rms(stroke_1dm(2,:));

figure()
plot(t,stroke_cascaded_standalone(2,t_start:t_end))
hold on;
plot(t,stroke_1dm(2,t_start:t_end),'LineWidth',1)
xlabel('Time (s)')
ylabel('Res. amp.')
title('Case 1')
legend('Case 1','Ref.','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/case_1.pdf -transparent

% export_fig ../plot/case_5.pdf -transparent

%% Case 3
rms_case3 = rms(stroke_parallel_standalone(2,:));
rms_ref = rms(stroke_1dm(2,:));

figure()
plot(t,stroke_parallel_standalone(2,t_start:t_end))
hold on;
plot(t,stroke_1dm(2,t_start:t_end),'LineWidth',1)
xlabel('Time (s)')
ylabel('Amp.')
title('Stroke')
legend('Case 3','Ref.','Interpreter','latex',Location='northwest');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/stroke_case_3.pdf -transparent

%% Case 4
rms_case3 = rms(stroke_parallel_integrated(2,:));
rms_ref = rms(stroke_1dm(2,:));

figure()
plot(t,stroke_parallel_integrated(2,t_start:t_end))
hold on;
plot(t,stroke_1dm(2,t_start:t_end),'LineWidth',1)
xlabel('Time (s)')
ylabel('Amp.')
title('Stroke')
legend('Case 3','Ref.','Interpreter','latex',Location='northwest');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
% export_fig ../plot/stroke_case_4.pdf -transparent
