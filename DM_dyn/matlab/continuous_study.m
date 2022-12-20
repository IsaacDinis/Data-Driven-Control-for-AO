clear;

Ts=1/1000;
s=tf('s');

%% Continuous time  model

w_DM = 20000*2*pi;

s_DM = 0.1;

DM = w_DM.^2/(s^2+2*s_DM*w_DM*s+w_DM.^2);

g = 0.5;
K = g/(Ts*s);
delay = exp(-s*650*1e-6);
%%
opts = bodeoptions('cstprefs');
opts.PhaseVisible = 'on';
opts.FreqUnits = 'Hz';
% opts.Title.String = 'Sensitivity function of the integrator with different gains';
opts.XLabel.Interpreter = 'latex';
opts.YLabel.Interpreter = 'latex';
opts.Title.Interpreter = 'latex';
opts.Title.FontSize = 15;
opts.YLabel.FontSize = 15;
opts.XLabel.FontSize = 15;
opts.TickLabel.FontSize = 15;


%% DM frequency response

%open loop
opts.Title.String = "DM dynamics Bode plot";
figure()
bodeplot(DM,opts);
grid on;
set(gcf, 'Position',  [100, 100, 700, 580])
% export_fig ../plot/DM_dyn.pdf -noinvert



%closed loop system low gain
g = 0.5;
K = g/(Ts*s);
opts.Title.String = "Closed-loop ($\mathcal{S}$) Bode plot, DM dynamics effect with low integrator gain";
opts.Ylim{1} = [-39 10];
opts.Ylim{2} = [-45 90];
opts.Xlim{1} = [1 1e4];
figure()
bodeplot(feedback(1,K*delay),feedback(1,K*DM*delay),opts);
grid on;
set(gcf, 'Position',  [100, 100, 700, 580])
legend('without DM dyn.','with DM dyn.','Interpreter','latex');
% export_fig ../plot/Closed_loop_low_gain.pdf -noinvert


%closed loop system high gain
g = 2;
K = g/(Ts*s);
opts.Title.String = "Closed-loop ($\mathcal{S}$) Bode plot, DM dynamics effect with high integrator gain";
opts.Ylim{1} = [-50 30];
opts.Ylim{2} = [-55 110];
opts.Xlim{1} = [1 1e4];
figure()
bodeplot(feedback(1,K*delay),feedback(1,K*DM*delay),opts);
grid on;
set(gcf, 'Position',  [100, 100, 700, 580])
legend('without DM dyn.','with DM dyn.','Interpreter','latex');
% export_fig ../plot/Closed_loop_high_gain.pdf -noinvert

%%
%open loop system low gain
g = 0.5;
K = g/(Ts*s);

delay = exp(-s*0*1e-6);
S1 = K*DM*delay;

delay = exp(-s*150*1e-6);
S2 = K*DM*delay;

delay = exp(-s*300*1e-6);
S3 = K*DM*delay;

delay = exp(-s*450*1e-6);
S4 = K*DM*delay;

delay = exp(-s*600*1e-6);
S5 = K*DM*delay;


figure()
opts.Ylim{1} = [-85 20];
opts.Ylim{2} = [-360 0];
opts.Title.String = "Open-loop Bode plot, stability margins depending on group delay";
bodeplot(S1,S2,S3,S4,S5,opts);
grid on;
set(gcf, 'Position',  [100, 100, 700, 580])
legend('0 $\mu$s','150 $\mu$s','300 $\mu$s','450 $\mu$s','600 $\mu$s','Interpreter','latex')

% export_fig ../plot/openloop_margin.pdf -noinvert

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

opts.Title.String = "Closed-loop ($\mathcal{S}$) Bode plot depending on group delay with high gain integrator";
opts.Xlim{1} = [10 5*1e3];
figure()
bodeplot(S1,S2,S3,S4,S5,S6,opts);

set(gcf, 'Position',  [100, 100, 700, 580])
legend('150 $\mu$s','250 $\mu$s','350 $\mu$s','450 $\mu$s','550 $\mu$s','650 $\mu$s', 'Interpreter','latex')
% export_fig ../plot/cl_delay_high.pdf -noinvert
%% low gain 

g = 0.5;
K = g/(Ts*s);

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

opts.Title.String = "Closed-loop ($\mathcal{S}$) Bode plot depending on group delay with low gain integrator";
opts.Xlim{1} = [10 5*1e3];

figure()
bodeplot(S1,S2,S3,S4,S5,opts);

set(gcf, 'Position',  [100, 100, 700, 580])
legend('150 $\mu$s','350 $\mu$s','650 $\mu$s','850 $\mu$s','1250 $\mu$s', 'Interpreter','latex')

export_fig ../plot/cl_delay_low.pdf -noinvert









