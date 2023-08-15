dist_matrix_path = '../data/single_mode_dist_nonoise_1_2000.mat';
dist = load(dist_matrix_path).data';
fs = 4000;
%% DM1
fs_DM1=1/1000;
z = tf('z',1/fs_DM1);
g = 0.5;
K_DM1 = d2d(g/(1-z^-1),1/fs);

%% DM2
z = tf('z',1/fs);
g = 0.5;
K_DM2 = g/(1-z^-1);
%%
[b_lp, a_lp] = butter(2,50/fs*2,'low');
figure()
freqz(b_lp,a_lp,[],fs)

[b_hp, a_hp] = butter(2,10/fs*2,'high');
figure()
freqz(b_hp,a_hp,[],fs)

figure()
freqz(b_lp.*b_hp,a_lp.*a_hp,[],fs)
%%
t = 0:1/fs:size(dist,1)/fs-1/fs;
sys = feedback(1,K_DM1*K_DM2);
y = lsim(sys,dist,t);                                                                          