s=tf('s');
Ts=1/10; %replace with desired value;
z = tf('z',Ts);
G_d = 1/(z-0.1);
G_c = d2c(G_d,'tustin');
G2_c = exp(-s*Ts)/(1-0.1*exp(-s*Ts));
step(G_d,G_c,G2_c)
bode(G_d,G_c,G2_c)