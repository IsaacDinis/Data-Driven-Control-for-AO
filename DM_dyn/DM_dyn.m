s=tf('s');
Ts=1/2000; %replace with desired value;
z = tf('z',Ts);

w_DM = 1200*2*pi;
s_dm = 0.1;
DM_c = w_DM.^2/(s^2+2*s_dm*w_DM*s+w_DM.^2);

DM_d = c2d(DM_c,Ts,'zoh');

figure()
bode(DM_d,DM_c)
legend()