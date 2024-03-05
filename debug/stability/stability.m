
fs = 1000; % sampling frequency
gain = 1; % integrator gain
RTC_delai = 0; % number of RTC delay frames 
leak = 0.01;
K = tf([gain,0],[1,-1+leak],1/fs); % integrator transfer function

dummy = zeros(1,1+RTC_delai);
dummy(1) = 1;
RTC = tf([1],dummy,1/fs); % RTC transfer function

WFS = tf([1,1],[2,0],1/fs); % WFS transfer function


% sys_ol = WFS*RTC*K; % open loop transfer function
sys_ol = RTC*K; % open loop transfer function
sys_cl = feedback(1,sys_ol); % closed loop transfer function

disp(['Eigenvalues closed-loop : ',...
    num2str(max(abs(eig(sys_cl)))), ' (stable if < 1)'])  

figure()
bodemag(sys_cl)
title('Closed loop rejection function')
grid on

% Open loop bode with margins
figure()
margin(sys_ol)


figure()
h = pzplot(sys_cl);
title('Closed loop zeros and poles')
% grid on

figure()
nyquist(sys_ol)
title('Open loop Nyquist')

figure()
step(sys_cl)