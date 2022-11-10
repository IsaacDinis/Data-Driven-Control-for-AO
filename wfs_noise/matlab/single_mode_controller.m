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

dist_matrix_path = '../data/single_mode_dist_noise.mat';
dist_matrix = load(dist_matrix_path).data';
dist_matrix = dist_matrix(2:end);
%%
fs = 1000;
bandwidth = 100;
order = 3;
max_control_gain = 20;
%%

window_size = 499;
n_average = 20;

[psd, f] = compute_psd(dist_matrix,n_average,window_size,fs);
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
W1 = W1/val;

W3 = tf(max_control_gain,1,1/fs);
Kdd =  ao_dd_controller(fs,w,order,W1,W3,'fusion');
Kdd_numerator = Kdd.Numerator{1};
Kdd_denominator = Kdd.Denominator{1};

toc

Kdd_matrix(:,1,:) = Kdd_numerator;
Kdd_matrix(:,2,:) = Kdd_denominator;

%% Simulation
g = 0.1;
K0 = tf([g,0],[1,-1],1/fs);

sys_dd = feedback(1,G*Kdd);
sys_int = feedback(1,G*K0);

t = 0:1/fs:size(dist_matrix,1)/fs-1/fs;

[ydd,~] = lsim(sys_dd,dist_matrix,t);
[yint,~] = lsim(sys_int,dist_matrix,t);

fprintf("rms datadriven = %f \nrms integrator = %f \n",rms(ydd(500:end)),rms(yint(500:end)));

%% Save controller coefficients
Kdd_matrix = reshape(Kdd_matrix,(max_order+1)*2,n_modes);
% % Kdd = circshift(Kdd,-2,2); % put tip tilt controllers at the end 
save('../Kdd_noise2','Kdd_matrix');
