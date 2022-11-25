%%
fs = 1000;
Ts = 1/fs;
s = tf('s');
z = tf('z',Ts);
Td = 650e-6; % AO loop delay 
g = 0.5;
%%
WFS_c = (1-exp(-s*Ts))/(s*Ts);
K_d = g/(1-z^-1);
RTC_c = d2c(K_d,'zoh')*exp(-s*Td);
K_c = d2c(K_d,'zoh'); 
K2_c = g/(1-exp(-s*Ts));
ZOH_c = (1-exp(-s*Ts))/(s*Ts);
%%
u = ones(1,length(t));
T = 1;
t = 0:1/fs:T-1/fs;


S1 = feedback(1,K_c*ZOH_c);
S2 = feedback(1,K2_c);
S3 = feedback(1,K_d);
figure()
step(S1,S2,S3);
figure()
lsim(S1,u,t)