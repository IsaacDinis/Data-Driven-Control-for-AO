% dist_matrix_path = '../data/single_mode_dist_nonoise_1_4000.mat';
% dist = load(dist_matrix_path).data';

%%
s=tf('s');
Ts=1/4000;
z = tf('z',Ts);


w_DM = 1200*2*pi;

s_dm = 0.01;
DM_c = w_DM.^2/(s^2+2*s_dm*w_DM*s+w_DM.^2);

Ts_K = 1/2000;
z_K = tf('z',Ts_K);
g = 0.5;
K = g/(1-z_K^-1);
K_sim = d2d(K,Ts,'tustin');

figure()
bode(K,K_sim);
legend();

delay = z_K^-1;

w_DM_d = 2/Ts*tan(w_DM*Ts/2);
DM_c_d = w_DM_d.^2/(s^2+2*s_dm*w_DM_d*s+w_DM_d.^2);
opt = c2dOptions('Method','tustin','PrewarpFrequency',w_DM);
DM_d_tustin = c2d(DM_c,Ts,opt);

DM_d_zoh = c2d(DM_c,Ts);
% figure()
% bode(DM_c,DM_d_tustin,DM_d_zoh);
% legend();
%%
t = 0:Ts:size(dist,1)*Ts-Ts;
% sys = feedback(1,K_sim*)