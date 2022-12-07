dist_matrix_path = '../data/single_mode_dist_nonoise_1_4000.mat';
dist = load(dist_matrix_path).data';

%%
s=tf('s');
Ts=1/4000;
z = tf('z',Ts);


w_DM = 1200*2*pi;

s_dm = 0.01;
DM_c = w_DM.^2/(s^2+2*s_dm*w_DM*s+w_DM.^2);


g = 0.1;
K = g/(1-z^-1);



delay = 1;

% w_DM_d = 2/Ts*tan(w_DM*Ts/2);
% DM_c_d = w_DM_d.^2/(s^2+2*s_dm*w_DM_d*s+w_DM_d.^2);
% opt = c2dOptions('Method','tustin','PrewarpFrequency',w_DM);
% DM_d_tustin = c2d(DM_c,Ts,opt);

DM_d_zoh = c2d(DM_c,Ts);
% figure()
% bode(DM_c,DM_d_tustin,DM_d_zoh);
% legend();

sys_ol = delay*K*DM_d_zoh;
sys_cl = feedback(1,sys_ol);
%%
figure()
% % bode(K2_c*ZOH_c,K_c,K_d)
h = bodeplot(sys_ol);
% setoptions(h,'FreqUnits','Hz','PhaseVisible','off');
setoptions(h,'FreqUnits','Hz');
% legend()
figure()
nyquist(sys_ol)

figure()
step(sys_cl)
%%
t = 0:Ts:size(dist,1)*Ts-Ts;

% sys = feedback(1,K_sim*)

