fs = 4000;
g = 0.5;
G = tf([1],[1,0,0],1/fs); % 2 samples delay
K = tf([g,0],[1,-1],1/fs);
sys = feedback(1,G*K);

figure()
bodemag(sys)