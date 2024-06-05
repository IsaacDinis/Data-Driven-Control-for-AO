fs = 4000;

s = tf('s');    % continuous frequency
z = tf('z',1/fs);
delay = z^-2; 
g = 0.5; % gain
K = g/(1-z^-1); % integrator
% K.InputDelay = 2;
% K.OutputDelay = 1;
K_c = d2c(K,'tustin')*exp(-2*1/fs*s);
T = 10e-3;



w = 4000*2*pi; % resonant frequency
zeta = 1;    % damping factor

DM_dyn = w.^2/(s^2+2*zeta*w*s+w.^2); % DM dynamics
% DM_dyn = s/s; % DM dynamics
sys_cl = feedback(1,DM_dyn*K_c);

fs_sim = 20000;
t = 0:1/fs_sim:T-1/fs_sim;

ol_step = step(DM_dyn,t);
cl_step = step(sys_cl,t);
cl_step_wo_dealay = step(feedback(1,K_c),t);

figure()
plot(t,ol_step);
hold on;
plot(t,cl_step);
plot(t,cl_step_wo_dealay);
legend('DM dynamics ol response','cl response with DM dynamics','cl response without DM dynamics')
xlabel('time (s)')
title('step response, continuous time')
make_it_nicer()

%%

DM_dyn_d = c2d(DM_dyn,1/fs,'foh'); % DM dynamics
% DM_dyn = s/s; % DM dynamics
sys_cl = feedback(1,DM_dyn_d*K*1/z*1/z);

fs_sim = 4000;
t = 0:1/fs_sim:T-1/fs_sim;

ol_step = step(DM_dyn_d,t);
cl_step = step(sys_cl,t);
cl_step_wo_dealay = step(feedback(1,K*1/z*1/z),t);

figure()
plot(t,ol_step);
hold on;
plot(t,cl_step);
plot(t,cl_step_wo_dealay);
legend('DM dynamics ol response','cl response with DM dynamics','cl response without DM dynamics')
xlabel('time (s)')
title('step response, discrete time')
make_it_nicer()