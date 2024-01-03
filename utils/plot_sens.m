fs = 1000;
K = tf([0.3],[1,-1],1/fs);
G = tf([1],[1,0,0],1/fs);
S = feedback(1,G*K);
max(abs(eig(S)))

K2 = tf([0.5,0.2],[1,-1],1/fs);

S2 = feedback(1,G*K2);
max(abs(eig(S2)))

figure()
bodemag(S,S2)
legend('integrateur simple', 'PI')
K3 = tf([0.5],[1,-0.99],1/fs);
% figure()
% bodemag(K,K3)

% figure()
% bodemag(feedback(G*K,G*K),feedback(G*K3,G*K3))

w = logspace(log10(0.00001),log10(pi*fs),10000);
[mag,phase,wout] = bode(S,w);
trapz(wout,log(mag))
[mag,phase,wout] = bode(S2,w);
trapz(wout,log(mag))

%%
s = tf('s');
plop = 1/(s*s);
% plop = 1/s;
K3 = c2d(plop,1/fs);
K3 = K3*0.0001*fs;
S3 = feedback(1,G*K3);
S4 = feedback(1,G*K*K*0.1);
max(abs(eig(S3)))
figure()
bodemag(S,S3)