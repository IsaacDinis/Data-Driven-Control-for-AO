%% load libraries
% clear all;
% tbxmanager restorepath
% addpath /home/isaac/mosek/9.3/toolbox/r2015a
% path_to_fusion = "/home/isaac/mosek/9.3/tools/platform/linux64x86/bin/mosek.jar";
path_to_fusion = "/home/isaac/mosek/10.1/tools/platform/linux64x86/bin/mosek.jar";
fusion_chk = contains(javaclasspath('-dynamic'), "mosek", 'IgnoreCase', true);
if strlength(path_to_fusion) && ~sum(fusion_chk)
    javaaddpath(path_to_fusion);
end 
%% Load data
mode_train = 2;
mode_test = 2;
RTC_delai = 2;
case_path = "../results/";
slopes_cl = fitsread(case_path+'int/saxoplus_KL_res.fits');
command_cl = fitsread(case_path+'int/saxoplus_KL_u.fits');

slopes_cl = slopes_cl(:,mode_train);
command_cl = command_cl(:,mode_train);
dist_matrix = command_cl(1:end-RTC_delai)+slopes_cl(1+RTC_delai:end);
dist_matrix = dist_matrix(100:end);
figure()
plot(dist_matrix)
%%
fs = 3000;
bandwidth = 180;
order = 50;
max_control_gain = 2;

%%

% window_size = 550*2*5;
% n_average = 20;

fft_size = 500;

% n_average = 15;
[psd,f] = compute_psd_welch(dist_matrix,fft_size,fs);

% psd(20:end) = psd(20);

% [psd, f] = periodogram(dist_matrix,[],length(dist_matrix),fs);
% [psd, f] =  pwelch(dist_matrix,400,200,400,fs);
% 
% psd =  fitsread('../python/psd.fits');
% psd = psd(:,mode);
% f =  fitsread('../python/freq.fits');
% f = f(2:end);
% psd = psd(2:end,:);

% psd(1:13) = psd(1:13)/8;
%%
psd(90:end) = mean(psd(90:end));
figure()
semilogx(f,10*log10(psd))
% psd(1:4) = psd(5);
G = tf([1],[1,0,0],1/fs); % 1 samples delay
% G = tf([1],[1,0],1/fs); % 1 samples delay

%% Load data
max_order =  max(order);
n_modes = 1;
Kdd_matrix = zeros(max_order+1,2,n_modes);

g = 0.256;
K0 = tf([g,0],[1,-1],1/fs);
S_int = feedback(1,G*K0);

tic
w = f*2*pi;
W1 = frd(psd,w,1/fs);
% W1 = frd(psd,w,1/fs)/abs(frd(S_int,w,'rad/s'));
val = freqresp(W1, bandwidth*2*pi);
% val = interp1(f,psd,bandwidth);
% val = 1;
W1 = W1/val;




% W13 = W1*W12;
W3 = tf(max_control_gain,1,1/fs);

% W32 = makeweight(0.001,200,10,1/fs);
% figure()
% bodemag(W1,W12,W13,W32)
% legend()

Kdd =  ao_dd_controller(fs,w,order,W1,W3,[],'fusion');
Kdd_numerator = Kdd.Numerator{1};
Kdd_denominator = Kdd.Denominator{1};

toc

% Kdd_matrix(:,1,:) = Kdd_numerator;
% Kdd_matrix(:,2,:) = Kdd_denominator;

%% Simulation

% dist_matrix_path = '../results/bright1_1_3ms/dcao/saxoplus_KL_u.fits';
% % dist_matrix_path = '../results/red3_05_55ms/dcao/saxoplus_KL_u.fits';
% % dist_matrix_path = '../data/single_mode_dist.mat';
% dist_matrix = fitsread(dist_matrix_path);
% dist_matrix = dist_matrix(:,mode);

K0 = tf([g,0],[1,-1],1/fs);

sys_dd = feedback(1,G*Kdd);
sys_int = feedback(1,G*K0);

t = 0:1/fs:size(dist_matrix,1)/fs-1/fs;

[ydd,~] = lsim(sys_dd,dist_matrix,t);
[yint,~] = lsim(sys_int,dist_matrix,t);

rms_dd = rms(ydd(500:end));
rms_int = rms(yint(500:end));

fprintf("rms datadriven = %f \nrms integrator = %f \n",rms_dd,rms_int);
gain = mean(100-rms_dd*100./rms_int);

[psd_dd, f] = compute_psd_fft(ydd,fft_size,fs);
[psd_int, f] = compute_psd_fft(yint,fft_size,fs);

figure()
semilogx(f,10*log10(psd_dd))
hold on;
semilogx(f,10*log10(psd_int))

legend('Data-driven','Integrator','Interpreter','latex','Location','southeast')
title('Residual PSD')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

grid()
make_it_nicer()

set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
% export_fig ../plot/residual_npnoise.pdf -transparent

%% Save controller coefficients
% Kdd_matrix = reshape(Kdd_matrix,(max_order+1)*2,n_modes);
% % Kdd = circshift(Kdd,-2,2); % put tip tilt controllers at the end 
% save('../Kdd_ProxCen','Kdd_matrix');
%%
S_dd = feedback(1,G*Kdd);
dummy = S_dd*W1*val;
figure()
bodemag(dummy,S_dd,W1)

S_int = feedback(1,G*K0);
dummy = S_int*W1*val;
figure()
bodemag(dummy,S_int,W1)

%%

Kdd_simp = tf(reduce(ss(Kdd),10));

S_dd_simp = feedback(1,G*Kdd_simp);

figure()
bodemag(S_dd_simp,S_dd);
figure()
pzmap(Kdd)

%%
[mag_dd,phase,wout] = bode(S_dd,w);
[mag_int,phase,wout] = bode(S_int,w);
mag_dd = 20*log10(squeeze(mag_dd));
mag_int = 20*log10(squeeze(mag_int));
f_bode = wout/2/pi;

[psd,f] = compute_psd_welch(dist_matrix,fft_size,fs);

figure()
semilogx(f,20*log10(val./psd),'color',[0.9290 0.6940 0.1250])
hold on;
semilogx(f_bode,mag_int,'color',[0 0.4470 0.7410])
semilogx(f_bode,mag_dd,'color',[0.8500 0.3250 0.0980])
xlabel('freq (Hz)')
ylabel('magnitude (dB)')
legend('inverse of dist. PSD','integrator','data-driven','Location','northwest')
title('Bode magnitude plot of sensitivity functions')
grid()
make_it_nicer()

set(gcf,'Position',[100 100 800 500])
set(gcf,'PaperType','A4')



figure()
semilogx(f,20*log10(psd_int))
hold on;
semilogx(f,20*log10(psd_dd))
xlabel('freq (Hz)')
ylabel('magnitude (dB)')
legend('integrator','data-driven','Location','northwest')
title('Tilt residual PSD')
grid()
make_it_nicer()

set(gcf,'Position',[100 100 800 500])
set(gcf,'PaperType','A4')

