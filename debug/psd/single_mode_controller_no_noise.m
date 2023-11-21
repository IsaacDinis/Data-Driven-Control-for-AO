%% load libraries
% clear all;
% tbxmanager restorepath
% addpath /home/isaac/mosek/9.3/toolbox/r2015a
path_to_fusion = "/home/isaac/mosek/10.0/tools/platform/linux64x86/bin/mosek.jar";
fusion_chk = contains(javaclasspath('-dynamic'), "mosek", 'IgnoreCase', true);
if strlength(path_to_fusion) && ~sum(fusion_chk)
    javaaddpath(path_to_fusion);
end 
%% Load data

dist_matrix_path = 'zernike_res.fits';
% dist_matrix_path = '../data/single_mode_dist.mat';
dist_matrix = fitsread(dist_matrix_path);
dist_matrix = dist_matrix(2:end,1);
figure()
plot(dist_matrix)
%%
fs = 2000;
bandwidth = 200;
order = 3;
max_control_gain = 0.8;

%%

% window_size = 550*2*5;
% n_average = 20;

window_size = 200;
n_average = 10;
nfft = 256;
[psd, f] = compute_psd2(dist_matrix,n_average,window_size,fs);
[psd_welch, f_welch] = pwelch(dist_matrix,400,200,400,fs);
% psd(1:13) = psd(1:13)/8;

% [psd, f] = periodogram(dist_matrix,[],[],fs);
figure()
semilogx(f,10*log10(psd))
hold on
semilogx(f_welch,10*log10(psd_welch))

% psd(1:4) = psd(5);
G = tf([1],[1,0,0],1/fs); % 1 samples delay

psd = psd_welch;
f = f_welch;
%% Load data
max_order =  max(order);
n_modes = 1;
Kdd_matrix = zeros(max_order+1,2,n_modes);


tic
w = f*2*pi;
W1 = frd(psd,w,1/fs);
val = freqresp(W1, bandwidth*2*pi);
% val = interp1(f,psd,bandwidth);
W1 = W1/val;




% W13 = W1*W12;
m = 0.1;
W3 = tf(max_control_gain,1,1/fs);
w32 = tf(1/m,1,1/fs);
% W32 = makeweight(0.001,200,10,1/fs);
% figure()
% bodemag(W1,W12,W13,W32)
% legend()

Kdd =  ao_dd_controller(fs,w,order,W1,W3,[],'fusion');
Kdd_numerator = Kdd.Numerator{1};
Kdd_denominator = Kdd.Denominator{1};

toc

Kdd_matrix(:,1,:) = Kdd_numerator;
Kdd_matrix(:,2,:) = Kdd_denominator;

%% Simulation
g = 0.5;
K0 = tf([g,0],[1,-1],1/fs);

sys_dd = feedback(1,Kdd*G);
sys_int = feedback(1,K0*G);

t = 0:1/fs:size(dist_matrix,1)/fs-1/fs;

[ydd,~] = lsim(sys_dd,dist_matrix,t);
[yint,~] = lsim(sys_int,dist_matrix,t);

rms_dd = rms(ydd(500:end));
rms_int = rms(yint(500:end));
rms_dist = rms(dist_matrix(500:end));

fprintf("rms disr = %f \n rms datadriven = %f \nrms integrator = %f \n",rms_dist,rms_dd,rms_int);
gain = mean(100-rms_dd*100./rms_int);

% [psd_dd, f] = compute_psd(ydd,n_average,window_size,fs);
% [psd_int, f] = compute_psd(yint,n_average,window_size,fs);

[psd_dd, f] = pwelch(ydd,400,200,400,fs);
[psd_int, f] = pwelch(yint,400,200,400,fs);
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
Kdd_matrix = reshape(Kdd_matrix,(max_order+1)*2,n_modes);
% % Kdd = circshift(Kdd,-2,2); % put tip tilt controllers at the end 
% save('../Kdd_ProxCen','Kdd_matrix');
%%
S_dd = feedback(1,G*Kdd);
dummy = S_dd*W1*val;
figure()
bodemag(dummy,S_dd,W1)