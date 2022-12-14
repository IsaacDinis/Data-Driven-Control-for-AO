clear;


Ts=1/4000;
s=tf('s');
z = tf('z',Ts);




%% bode options
opts = bodeoptions('cstprefs');
opts.PhaseVisible = 'on';
opts.FreqUnits = 'Hz';
% opts.Title.String = 'Sensitivity function of the integrator with different gains';
opts.XLabel.Interpreter = 'latex';
opts.YLabel.Interpreter = 'latex';
opts.Title.Interpreter = 'latex';
opts.Title.FontSize = 15;
opts.YLabel.FontSize = 15;
opts.XLabel.FontSize = 15;
opts.TickLabel.FontSize = 15;
%% Continuous time  model

w_DM = 1200*2*pi;

s_DM = 0.1;

DM_c = w_DM.^2/(s^2+2*s_DM*w_DM*s+w_DM.^2);

K_DM_c = (s^2+2*s_DM*w_DM*s+w_DM.^2)/(s^2+2*w_DM*s+w_DM.^2);



K_int_c = 0.5/(1e-3*s);
delay = exp(-s*2*Ts);



%% Continuous time plots

opts.YlimMode = 'auto';
opts.XlimMode = 'auto';

opts.Title.String = "Bode plot of DM damper";
figure()
bode(DM_c,K_DM_c,DM_c*K_DM_c,opts);
set(gcf, 'Position',  [100, 100, 700, 580])
legend('DM','damper','DM$\times$damper','Interpreter','latex');
grid on;
% export_fig ../plot/Damper_ol.pdf -noinvert

opts.Title.String = "Closed-loop ($\mathcal{S}$) Bode plot";
figure()
bode(feedback(1,DM_c),feedback(1,DM_c*K_DM_c),opts);
set(gcf, 'Position',  [100, 100, 700, 580])
legend('w/o damper','w/ damper','Interpreter','latex');
grid on;
% export_fig ../plot/damper_cl.pdf -noinvert

[y1,tOut] = step(feedback(1,DM_c));
[y2] = step(feedback(1,DM_c*K_DM_c),tOut);
figure()
plot(tOut,y1)
hold on;
plot(tOut,y2)
xlim([tOut(1),tOut(end)]);
title("Step response")
xlabel("Time (s)")
ylabel("Amplitude")
legend('w/o damper','w/ damper','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 290])
% export_fig ../plot/damper_step.pdf -noinvert

%% Same but  integrator

opts.Title.String = "Open-loop Bode plot with integrator";
figure()
bode(DM_c*K_int_c,DM_c*K_DM_c*K_int_c,opts);
set(gcf, 'Position',  [100, 100, 700, 580])
legend('w/o damper','w/ damper','Interpreter','latex');
grid on;
% export_fig ../plot/damper_ol_int.pdf -noinvert

opts.Title.String = "Closed-loop ($\mathcal{S}$) Bode plot with integrator";
figure()
bode(feedback(1,DM_c*K_int_c),feedback(1,DM_c*K_DM_c*K_int_c),opts);
set(gcf, 'Position',  [100, 100, 700, 580])
legend('w/o damper','w/ damper','Interpreter','latex');
grid on;
% export_fig ../plot/damper_cl_int.pdf -noinvert

[y1,tOut] = step(feedback(1,DM_c*K_int_c));
[y2] = step(feedback(1,DM_c*K_DM_c*K_int_c),tOut);
figure()
plot(tOut,y1)
hold on;
plot(tOut,y2)
xlim([tOut(1),tOut(end-50)]);
title("Step response with integrator")
xlabel("Time (s)")
ylabel("Amplitude")
legend('w/o damper','w/ damper','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 290])
% export_fig ../plot/damper_step_int.pdf -noinvert
%% Same but with delay and integrator

sys_ol = DM_c*K_DM_c;

% sys_cl = feedback(1,sys_ol);

% figure()
% bode(DM_c,K_DM_c,DM_c*K_DM_c);
% 
% figure()
% bode(sys_ol*K_int_c);

opts.Ylim{1} = [-100 20];
opts.Ylim{2} = [-720 0];
opts.Xlim{1} = [10 1e4];
opts.Title.String = "Open-loop Bode plot with integrator and group delay";
figure()
bode(DM_c*K_int_c*delay,DM_c*K_DM_c*K_int_c*delay,opts);
set(gcf, 'Position',  [100, 100, 700, 580])
legend('w/o damper','w/ damper','Interpreter','latex');
grid on;
% export_fig ../plot/damper_ol_int_delay.pdf -noinvert

opts.YlimMode = 'auto';
opts.XlimMode = 'auto';

opts.Title.String = "Closed-loop ($\mathcal{S}$) Bode plot with integrator and group delay";
figure()
bode(feedback(1,DM_c*delay*K_int_c),feedback(1,DM_c*K_DM_c*delay*K_int_c),opts);
set(gcf, 'Position',  [100, 100, 700, 580])
legend('w/o damper','w/ damper','Interpreter','latex');
grid on;

% export_fig ../plot/damper_cl_int_delay.pdf -noinvert

tOut = 0:1/4000:6e-3;
[y1] = step(feedback(1,DM_c*K_int_c*delay),tOut);
[y2] = step(feedback(1,DM_c*K_DM_c*K_int_c*delay),tOut);
figure()
plot(tOut,y1)
hold on;
plot(tOut,y2)
title("Step response with integrator and group delay")
xlabel("Time (s)")
ylabel("Amplitude")
legend('w/o damper','w/ damper','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 290])
% export_fig ../plot/damper_cl_int_delay.pdf -noinvert
%% Discrete time components

DM_d = c2d(DM_c,Ts);

opt = c2dOptions('Method','tustin','PrewarpFrequency',w_DM);

K_DM_d = c2d(K_DM_c,Ts,opt);


%% Case 1 no DM dynamics no photon noise
K_int_d = 0.5/(1-z^-1);


% figure()
% bode(sys_d_ol);
% 
% figure()
% step(sys_d_cl);
% 
opts.Xlim{1} = [1 2e3];
opts.YlimMode = 'auto';
% opts.XlimMode = 'auto';
opts.Title.String = "Closed-loop ($\mathcal{S}$) Bode plot with high gain integrator";
figure()
bode(feedback(1,K_int_d*z^-2),feedback(1,K_int_d*DM_d*z^-2),feedback(1,K_int_d*z^-2*K_DM_d),feedback(1,Kdd*z^-2*K_DM_d),opts);
set(gcf, 'Position',  [100, 100, 700, 580])
legend('w/o DM dyn.','w/o damper','w/ damper','data-driven','Interpreter','latex');
grid on;

% export_fig ../plot/damper_cl_int_delay.pdf -noinvert

dist_matrix_path = '../data/single_mode_dist_nonoise_1_4000.mat';
dist = load(dist_matrix_path).data';
t = 0:Ts:size(dist,1)*Ts-Ts;

y1 = lsim(feedback(1,K_int_d*z^-2),dist,t);
y2 = lsim(feedback(1,K_int_d*DM_d*z^-2),dist,t);
y3 = lsim(feedback(1,K_int_d*z^-2*K_DM_d),dist,t);

y4 = lsim(feedback(1,Kdd*K_DM_d*z^-2),dist,t);

rms_y1 = rms(y1);
rms_y2 = rms(y2);
rms_y3 = rms(y3);
rms_y4 = rms(y4);
% figure()
% bode(c2d(DM_c*delay,Ts),c2d(DM_c,Ts)*z^-1)


%%

K_int_d = 0.1/(1-z^-1);


% figure()
% bode(sys_d_ol);
% 
% figure()
% step(sys_d_cl);
% 
opts.Xlim{1} = [1 2e3];
opts.YlimMode = 'auto';
% opts.XlimMode = 'auto';
opts.Title.String = "Closed-loop ($\mathcal{S}$) Bode plot with low gain integrator";
figure()
bode(feedback(1,K_int_d*z^-2),feedback(1,K_int_d*DM_d*z^-2),feedback(1,K_int_d*z^-2*K_DM_d),opts);
set(gcf, 'Position',  [100, 100, 700, 580])
legend('w/o DM dyn.','w/o damper','w/ damper','Interpreter','latex');
grid on;
% 
% figure()
% pzplot(sys_d_cl)

dist_matrix_path = '../data/single_mode_dist_ProxCen_1_4000.mat';
dist = load(dist_matrix_path).data';
t = 0:Ts:size(dist,1)*Ts-Ts;

y1 = lsim(feedback(1,K_int_d*z^-2),dist,t);
y2 = lsim(feedback(1,K_int_d*DM_d*z^-2),dist,t);
y3 = lsim(feedback(1,K_int_d*z^-2*K_DM_d),dist,t);
rms_y1 = rms(y1);
rms_y2 = rms(y2);
rms_y3 = rms(y3);