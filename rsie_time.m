Ts=1/4000; % simulation sampling frequency
Ts_K = 1/500; % woofer sampling frequency​
z = tf('z',Ts); % simulation frequency
s = tf('s');    % continuous frequency
z_K = tf('z',Ts_K); % controller frequency​
w = 1200*2*pi; % resonant frequency
zeta = 0.1;    % damping factor​
DM_dyn_c = w.^2/(s^2+2*zeta*w*s+w.^2); % DM dynamics​
g = 0.5; % gain
K = g/(1-z_K^-1); % integrator
delay = z_K^-2;​
% convert to simulation space
K_sim = d2d(K,Ts,'tustin');
DM_dyn_d = c2d(DM_dyn_c,Ts,'tustin');
delay_sim = z^-(2*ceil(Ts_K/Ts));​
% Sensitivity functions
S_K = feedback(1,K*delay);
S_sim_no_dyn = feedback(1,K_sim*delay_sim);
S_sim = feedback(1,K_sim*delay_sim*DM_dyn_d);​
% plot step response of DM dynamics
figure()
step(DM_dyn_d)
hold on;
step(DM_dyn_c)
legend('continuous','discrete')
title('step response of DM dynamics')​
% plot closed loop step response 
figure()
step(S_K)
hold on;
step(S_sim_no_dyn)
step(S_sim)
legend('controller framerate','higher framerate without DM dynamics','higher framerate with DM dynamics')
title('closed loop step response')









