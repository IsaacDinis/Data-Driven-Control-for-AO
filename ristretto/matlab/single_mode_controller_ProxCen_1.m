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

dist_matrix_path = '../data/single_mode_dist_ProxCen_1.mat';
dist_matrix = load(dist_matrix_path).data';
dist_matrix = dist_matrix(2:end);
%%
fs = 1000;
bandwidth = 400;
order = 8;
max_control_gain = 0.1;
noise_red = 4;
%%

window_size = 49*5;
n_average = 5*2*20;

[psd, f] = compute_psd(dist_matrix,n_average,window_size,fs);
[psd, f] = periodogram(dist_matrix,[],[],fs);
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

% W12 =  (makeweight(1900,35*2*pi*2,0.9,1/fs));
% W12 =  (makeweight(0.025,35*2*pi*2,1.1,1/fs));

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

Kdd_matrix(:,1,:) = Kdd_numerator;
Kdd_matrix(:,2,:) = Kdd_denominator;

%% Simulation
g = 0.3;
K0 = tf([g,0],[1,-1],1/fs);

sys_dd = feedback(1,Kdd*G);
sys_int = feedback(1,K0*G);

t = 0:1/fs:size(dist_matrix,1)/fs-1/fs;

[ydd,~] = lsim(sys_dd,dist_matrix,t);
[yint,~] = lsim(sys_int,dist_matrix,t);

fprintf("rms datadriven = %f \nrms integrator = %f \n",rms(ydd(500:end)),rms(yint(500:end)));

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