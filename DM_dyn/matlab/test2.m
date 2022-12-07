clear;
Ts=1/5000;
s=tf('s');


w_DM = 1200*2*pi;

s_dm = 0.01;

DM_c = w_DM.^2/(s^2+2*s_dm*w_DM*s+w_DM.^2);

K = (w_DM.^2/(s^2+2*w_DM*s+w_DM.^2))/DM_c;
K2 = 1/s;

delay = exp(-s*Ts);

wc = 0.5; 
Gd = wc/s;

% K = loopsyn(DM_c,Gd);
sys_ol = DM_c*K*K2*delay;

sys_cl = feedback(1,sys_ol);

figure()
bode(sys_ol);

figure()
step(sys_cl);

figure()
bode(sys_cl);

%%
Ts=1/5000;
z = tf('z',Ts);

DM_d = c2d(DM_c*delay,Ts);
opt = c2dOptions('Method','tustin','PrewarpFrequency',w_DM);
K_d_dm = c2d(K,Ts,opt);
K_d_int = 0.1/(1-z^-1);
sys_d_ol = DM_d*K_d_dm*K_d_int;

sys_d_cl = feedback(1,sys_d_ol);

figure()
bode(sys_d_ol);

figure()
step(sys_d_cl);

figure()
bode(sys_d_cl);

% figure()
% bode(c2d(DM_c*delay,Ts),c2d(DM_c,Ts)*z^-1)