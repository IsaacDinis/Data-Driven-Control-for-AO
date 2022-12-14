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
ZOH_c = (1-exp(-s*Ts))/(Ts*s);

w_DM = 1200*2*pi;
s_dm = 0.1;
DM_c = w_DM.^2/(s^2+2*s_dm*w_DM*s+w_DM.^2);
DM_d = c2d(DM_c,Ts,'tustin');
%%
T = 0.02;
fs_sim = 4*fs; 
t = 0:1/fs_sim:T-1/fs_sim;

u = ones(1,length(t));

S1 = feedback(1,K_c);
S2 = feedback(1,K_c*DM_c);
S3 = feedback(1,K_d*DM_d);

% figure()
% step(S1,S2,S3);
% legend()

figure()
bodemag(S1,S2,S3)

figure()
bodemag(DM_c)
legend()

% 
% figure()
% lsim(S1,u,t)
% https://eng.libretexts.org/Bookshelves/Industrial_and_Systems_Engineering/Book%3A_Introduction_to_Control_Systems_(Iqbal)/07%3A_Design_of_Sampled-Data_Systems/7.02%3A_Pulse_Transfer_Function