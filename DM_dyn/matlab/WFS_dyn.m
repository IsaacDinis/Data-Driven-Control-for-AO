s=tf('s');
Ts=1/2000; %replace with desired value;
z = tf('z',Ts);

WFS_c = (1-exp(-Ts*s))/(Ts*s);
WFS_d = (z+1)/(2*z);
WFS2_d = c2d(WFS_c,Ts,'tustin');

figure()
bode(WFS_c,WFS_d,WFS2_d)
legend()

figure()
impulse(WFS_c,WFS_d,WFS2_d)
legend()