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

dist_matrix_path = '../data/single_mode_dist_PDS70_1.mat';
dist_matrix = load(dist_matrix_path).data';
dist_matrix = dist_matrix(2:end);
%%
fs = 1000;
bandwidth = 50;
order = 5;
max_control_gain = 0.1;
noise_red = 4;
%%

window_size = 499*4;
n_average = 5*5*2;

[psd, f] = compute_psd(dist_matrix,n_average,window_size,fs);
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
W12 =  (makeweight(0.025,35*2*pi*2,4,1/fs));

W13 = W1*W12;
W3 = tf(max_control_gain,1,1/fs);

W32 = makeweight(0.001,200,50,1/fs);
figure()
bodemag(W1,W12,W13,W32)
legend()

Kdd =  ao_dd_controller(fs,w,order,W1,W3,W32,'fusion');
Kdd_numerator = Kdd.Numerator{1};
Kdd_denominator = Kdd.Denominator{1};

toc

Kdd_matrix(:,1,:) = Kdd_numerator;
Kdd_matrix(:,2,:) = Kdd_denominator;

%% Simulation
g = 0.3;
K0 = tf([g,0],[1,-1],1/fs);

sys_dd = feedback(1,G*Kdd);
sys_int = feedback(1,G*K0);

t = 0:1/fs:size(dist_matrix,1)/fs-1/fs;

[ydd,~] = lsim(sys_dd,dist_matrix,t);
[yint,~] = lsim(sys_int,dist_matrix,t);

fprintf("rms datadriven = %f \nrms integrator = %f \n",rms(ydd(500:end)),rms(yint(500:end)));

[psd_dd, f] = compute_psd(ydd,n_average,window_size,fs);
[psd_int, f] = compute_psd(yint,n_average,window_size,fs);

figure()
semilogx(f,10*log10(psd_dd))
hold on;
semilogx(f,10*log10(psd_int))

legend('dd','int')
xlabel('frequency [Hz]')
ylabel('magnitude [dB]')

figure()
plot(f,10*log10(psd_dd))
hold on;
semilogx(f,10*log10(psd_int))

plot('dd','int')
xlabel('frequency [Hz]')
ylabel('magnitude [dB]')

%% Save controller coefficients
% Kdd_matrix = reshape(Kdd_matrix,(max_order+1)*2,n_modes);
% % % Kdd = circshift(Kdd,-2,2); % put tip tilt controllers at the end 
% save('../Kdd_noise2','Kdd_matrix');
