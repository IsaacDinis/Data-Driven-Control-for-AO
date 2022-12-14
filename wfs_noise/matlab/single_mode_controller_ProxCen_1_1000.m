%% load libraries
% clear all;
% tbxmanager restorepath
% addpath /home/isaac/mosek/9.3/toolbox/r2015a
path_to_fusion = "/home/isaac/mosek/9.3/tools/platform/linux64x86/bin/mosek.jar";
fusion_chk = contains(javaclasspath('-dynamic'), "mosek", 'IgnoreCase', true);
if strlength(path_to_fusion) && ~sum(fusion_chk)
    javaaddpath(path_to_fusion);
end 
%% Load data

% dist_matrix_path = '../data/single_mode_dist_ProxCen_1.mat';
dist_matrix_path = '../data/single_mode_dist.mat';
dist_matrix = load(dist_matrix_path).data';
dist_matrix = dist_matrix(2:end);
%%
fs = 2000;
bandwidth = 400;
order = 2;
max_control_gain = 0.1;
noise_red = 4;
%%

% window_size = 550*2*5;
% n_average = 20;

window_size = 999;
n_average = 20;

[psd, f] = compute_psd(dist_matrix,n_average,window_size,fs);
% [psd, f] = periodogram(dist_matrix,[],[],fs);
figure()
semilogx(f,10*log10(psd))
% psd(1:4) = psd(5);
G = tf([1],[1,0,0],1/fs); % 2 samples delay


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

W12 = frd(makeweight(10^(54.5/20),[152 10^(34/20)],10^(20/20),1/fs,2),w);

W4 = makeweight(10^(0/20),[50,10^(1/20)],10^(30/20),1/fs,2);
plop = W1*W4/freqresp(W1*W4, bandwidth*2*pi);

val = freqresp(W12, bandwidth*2*pi);

W12 = W12/val;

figure()
bodemag(W1,W12,W4,plop)

W3 = tf(max_control_gain,1,1/fs);

Kdd =  ao_dd_controller(fs,w,order,plop,W3,W4,'fusion');
Kdd_numerator = Kdd.Numerator{1};
Kdd_denominator = Kdd.Denominator{1};

toc

Kdd_matrix(:,1,:) = Kdd_numerator;
Kdd_matrix(:,2,:) = Kdd_denominator;

%% Simulation
g = 0.2;
K0 = tf([g,0],[1,-1],1/fs);

sys_dd = feedback(1,Kdd*G);
sys_int = feedback(1,K0*G);

t = 0:1/fs:size(dist_matrix,1)/fs-1/fs;

[ydd,~] = lsim(sys_dd,dist_matrix,t);
[yint,~] = lsim(sys_int,dist_matrix,t);

rms_dd = rms(ydd(500:end));
rms_int = rms(yint(500:end));

fprintf("rms datadriven = %f \nrms integrator = %f \n",rms_dd,rms_int);
gain = mean(100-rms_dd*100./rms_int);

[psd_dd, f] = compute_psd(ydd,n_average,window_size,fs);
[psd_int, f] = compute_psd(yint,n_average,window_size,fs);

figure()
semilogx(w,10*log10(psd_dd))
hold on;
semilogx(w,10*log10(psd_int))

legend('dd','int')
xlabel('frequency [Hz]')
ylabel('magnitude [dB]')

%% Save controller coefficients
Kdd_matrix = reshape(Kdd_matrix,(max_order+1)*2,n_modes);
% % Kdd = circshift(Kdd,-2,2); % put tip tilt controllers at the end 
save('../Kdd_ProxCen','Kdd_matrix');
%%
S_dd = feedback(1,G*Kdd);
dummy = S_dd*W1*val;
figure()
bodemag(dummy,S_dd,W1)