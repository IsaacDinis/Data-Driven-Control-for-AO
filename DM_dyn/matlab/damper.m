clear;

dist_matrix_path = '../data/single_mode_dist_ProxCen_1_4000.mat';
dist = load(dist_matrix_path).data';

Ts=1/4000;
s=tf('s');
z = tf('z',Ts);

t = 0:Ts:size(dist,1)*Ts-Ts;


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

K_int_c = 0.1/(Ts*s);

K_int_c2 = 0.5/(1e-3*s);
delay = exp(-s*2*Ts);



%% Continuous time plots



opts.Title.String = "Bode plot of DM damper";
figure()
bode(DM_c,K_DM_c,DM_c*K_DM_c,opts);
set(gcf, 'Position',  [100, 100, 700, 580])
legend('DM','damper','DM$\times$damper','Interpreter','latex');
grid on;

opts.Title.String = "Closed-loop Bode plot";
figure()
bode(feedback(1,DM_c),feedback(1,DM_c*K_DM_c),opts);
set(gcf, 'Position',  [100, 100, 700, 580])
legend('without damper','with damper','Interpreter','latex');
grid on;

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
legend('without damper','with damper','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 290])

%% Same but  integrator


opts.Title.String = "Closed-loop Bode plot with integrator";
figure()
bode(feedback(1,DM_c*K_int_c),feedback(1,DM_c*K_DM_c*K_int_c),opts);
set(gcf, 'Position',  [100, 100, 700, 580])
legend('without damper','with damper','Interpreter','latex');
grid on;

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
legend('without damper','with damper','Interpreter','latex');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 290])
%% Same but with delay and integrator

sys_ol = DM_c*K_DM_c;

% sys_cl = feedback(1,sys_ol);

% figure()
% bode(DM_c,K_DM_c,DM_c*K_DM_c);
% 
% figure()
% bode(sys_ol*K_int_c);

figure()
bode(feedback(1,DM_c*delay*K_int_c),feedback(1,DM_c*K_DM_c*delay*K_int_c));

figure()
step(feedback(1,DM_c*delay*K_int_c),feedback(1,DM_c*K_DM_c*delay*K_int_c));

%% Discrete time components

DM_d = c2d(DM_c,Ts);

opt = c2dOptions('Method','tustin','PrewarpFrequency',w_DM);

K_DM_d = c2d(K_DM_c,Ts,opt);


%% Case 1 no DM dynamics
K_int_d = 0.2/(1-z^-1);

sys_d_ol = K_int_d*z^-2;

sys_d_cl = feedback(1,sys_d_ol);

% figure()
% bode(sys_d_ol);
% 
% figure()
% step(sys_d_cl);
% 
% figure()
% bode(sys_d_cl);
% 
% figure()
% pzplot(sys_d_cl)

y = lsim(sys_d_cl,dist,t);
rms(y)
% figure()
% bode(c2d(DM_c*delay,Ts),c2d(DM_c,Ts)*z^-1)

%% Case 2 DM dynamics without damping

K_int_d = 0.2/(1-z^-1);

sys_d_ol = DM_d*K_int_d*z^-2;

sys_d_cl = feedback(1,sys_d_ol);

figure()
bode(sys_d_ol);

figure()
step(sys_d_cl);

figure()
bode(sys_d_cl);

figure()
pzplot(sys_d_cl)

y = lsim(sys_d_cl,dist,t);
rms(y)

%%  Case 3 DM dynamics with damping

K_int_d = 0.2/(1-z^-1);

sys_d_ol = DM_d*K_DM_d*K_int_d*z^-2;

sys_d_cl = feedback(1,sys_d_ol);

figure()
bode(sys_d_ol);

figure()
step(sys_d_cl);

figure()
bode(sys_d_cl);

figure()
pzplot(sys_d_cl)

y = lsim(sys_d_cl,dist,t);
rms(y)


%% Case 4 DM dynamics with data-driven controller