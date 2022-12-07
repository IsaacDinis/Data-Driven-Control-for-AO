%%
rms_dm_delay_2000 = 4.7;
rms_nodm_delay_2000 = 2.7;

rms_dm_nodelay_2000 = 3.6;
rms_nodm_nodelay_2000 = 1.5;

rms_dm_delay_4000 = 2.5;
rms_nodm_delay_4000 = 0.98;

rms_dm_nodelay_4000 = "unstable";
rms_nodm_nodelay_4000 = 0.25;
%%
s=tf('s');
Ts=1/2000;
z = tf('z',Ts);

sys = (1-exp(-s*Ts))/s;
ZOH_c = (1-exp(-s*Ts))/(Ts*s);

g = 0.5;
K2_c = g/(1-exp(-s*Ts));
K_d = g/(1-z^-1);
K_c = d2c(K_d,'zoh');

w_DM = 1200*2*pi;
s_dm = 0.1;
DM_c = w_DM.^2/(s^2+2*s_dm*w_DM*s+w_DM.^2);

delay = exp(-s*650*1e-6);


sys = ZOH_c*DM_c*delay*g/(s*Ts);
sys_no_delay = ZOH_c*DM_c*g/(s*Ts);
sys_no_DM = ZOH_c*delay*g/(s*Ts);
sys_no_DM_no_delay

% S1 = feedback(1,K2_c*DM_c*plop);
S1 = feedback(1,sys);
S2 = feedback(1,sys_no_delay);
S3 = feedback(1,sys_no_DM);
S4 = feedback(1,sys_no_DM_no_delay);

% S3 = feedback(1,K2_c);
figure()
% % bode(K2_c*ZOH_c,K_c,K_d)
h = bodeplot(S1);
setoptions(h,'FreqUnits','Hz','PhaseVisible','off');
% legend()
%%
t_sim = 0:1/10000:1;
figure()
step(S1,t_sim)
legend()

figure()
h = bodeplot(DM_c);
setoptions(h,'FreqUnits','Hz');



%%
figure()
h = bodeplot(S1,S2);
setoptions(h,'FreqUnits','Hz','PhaseVisible','off');
legend('650 us delay', 'no delay')

figure()
h = bodeplot(S1,S3);
setoptions(h,'FreqUnits','Hz','PhaseVisible','off');
legend('with DM dynamics', 'without DM dynamics')