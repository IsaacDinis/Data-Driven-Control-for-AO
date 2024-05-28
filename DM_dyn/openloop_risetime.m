s = tf('s');    % continuous frequency



w = 4000*2*pi; % resonant frequency
zeta = 0.99;    % damping factor

DM_dyn_c = w.^2/(s^2+2*zeta*w*s+w.^2); % DM dynamics
figure()
step(DM_dyn_c)