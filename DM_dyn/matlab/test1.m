clear;


Ts=1/1000;
s=tf('s');
z = tf('z',Ts);

%%
K_int_d = 0.5/(1-z^-1);






dist_matrix_path = '../data/single_mode_dist_nonoise_1_4000.mat';
dist = load(dist_matrix_path).data';
dist = dist(1:4:end);
t = 0:Ts:size(dist,1)*Ts-Ts;
t = 0:Ts:0.01;
y1 = lsim(feedback(1,K_int_d*z^-2),dist(1:length(t)),t);

rms_y1 = rms(y1);

K_int_d = 0.2/(1-z^-1);
dist_matrix_path = '../data/single_mode_dist_ProxCen_1_4000.mat';
dist = load(dist_matrix_path).data';
dist = dist(1:4:end);
t = 0:Ts:size(dist,1)*Ts-Ts;
t = 0:Ts:0.04;

y2 = lsim(feedback(1,K_int_d*z^-2),dist(1:length(t)),t);

rms_y2 = rms(y2);
