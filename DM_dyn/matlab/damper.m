clear;
Ts=1/4000;
s=tf('s');
z = tf('z',Ts);

%% Continuous time  model

w_DM = 1200*2*pi;

s_DM = 0.1;

DM_c = w_DM.^2/(s^2+2*s_DM*w_DM*s+w_DM.^2);
K_DM_c = (s^2+2*s_DM*w_DM*s+w_DM.^2)/(s^2+2*w_DM*s+w_DM.^2);

K_int_c = 0.1/(Ts*s);

delay = exp(-s*2*Ts);



%% Continuous time plots

sys_ol = DM_c*K_DM_c;

% sys_cl = feedback(1,sys_ol);

figure()
bode(DM_c,K_DM_c,DM_c*K_DM_c);
% 
% figure()
% bode(sys_ol*K_int_c);

figure()
bode(feedback(1,DM_c),feedback(1,DM_c*K_DM_c));

figure()
step(feedback(1,DM_c),feedback(1,DM_c*K_DM_c));
%% Same but  integrator

sys_ol = DM_c*K_DM_c;

% sys_cl = feedback(1,sys_ol);

% figure()
% bode(DM_c,K_DM_c,DM_c*K_DM_c);
% 
figure()
bode(DM_c*K_int_c);

figure()
bode(feedback(1,DM_c*K_int_c),feedback(1,DM_c*K_DM_c*K_int_c));

figure()
step(feedback(1,DM_c*K_int_c),feedback(1,DM_c*K_DM_c*K_int_c));
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
K_int_d = 0.5/(1-z^-1);

sys_d_ol = K_int_d*z^-2;

sys_d_cl = feedback(1,sys_d_ol);

figure()
bode(sys_d_ol);

figure()
step(sys_d_cl);

figure()
bode(sys_d_cl);

figure()
pzplot(sys_d_cl)

% figure()
% bode(c2d(DM_c*delay,Ts),c2d(DM_c,Ts)*z^-1)

%% Case 2 DM dynamics without damping

K_int_d = 0.5/(1-z^-1);

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

%%  Case 3 DM dynamics with damping

K_int_d = 0.5/(1-z^-1);

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

%% Case 4 DM dynamics with data-driven controller
