clear;

Ts=1/1000;
s=tf('s');

%% Continuous time  model

w_DM = 1200*2*pi;

s_DM = 0.1;

DM = w_DM.^2/(s^2+2*s_DM*w_DM*s+w_DM.^2);

g = 0.5;
K = g/(Ts*s);
delay = exp(-s*650*1e-6);

%% DM frequency response

%open loop
figure()
h = bodeplot(DM);
setoptions(h,'FreqUnits','Hz');

%closed loop system low gain
g = 0.5;
K = g/(Ts*s);
figure()
h = bodeplot(feedback(1,K*DM*delay),feedback(1,K*delay));

%closed loop system high gain
g = 2;
K = g/(Ts*s);
figure()
h = bodeplot(feedback(1,K*DM*delay),feedback(1,K*delay));
setoptions(h,'FreqUnits','Hz');

%% impact of delay 
%high gain
g = 2;
K = g/(Ts*s);

delay = exp(-s*150*1e-6);
S1 = feedback(1,K*DM*delay);

delay = exp(-s*250*1e-6);
S2 = feedback(1,K*DM*delay);

delay = exp(-s*350*1e-6);
S3 = feedback(1,K*DM*delay);

delay = exp(-s*450*1e-6);
S4 = feedback(1,K*DM*delay);

delay = exp(-s*550*1e-6);
S5 = feedback(1,K*DM*delay);

delay = exp(-s*650*1e-6);
S6 = feedback(1,K*DM*delay);

h = bodeplot(S1,S2,S3,S4,S5,S6);
setoptions(h,'FreqUnits','Hz');
% low gain 

g = 0.5;
K = g/(Ts*s);
figure()

delay = exp(-s*150*1e-6);
S1 = feedback(1,K*DM*delay);

delay = exp(-s*350*1e-6);
S2 = feedback(1,K*DM*delay);

delay = exp(-s*650*1e-6);
S3 = feedback(1,K*DM*delay);

delay = exp(-s*850*1e-6);
S4 = feedback(1,K*DM*delay);

delay = exp(-s*1250*1e-6);
S5 = feedback(1,K*DM*delay);

figure()
h = bodeplot(S1,S2,S3,S4,S5);
setoptions(h,'FreqUnits','Hz');










