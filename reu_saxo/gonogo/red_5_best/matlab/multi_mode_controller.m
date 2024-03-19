%% load libraries
path_to_fusion = "/home/isaac/mosek/9.3/tools/platform/linux64x86/bin/mosek.jar";
% path_to_fusion = "/home/isaac/mosek/10.0/tools/platform/linux64x86/bin/mosek.jar";
fusion_chk = contains(javaclasspath('-dynamic'), "mosek", 'IgnoreCase', true);
if strlength(path_to_fusion) && ~sum(fusion_chk)
    javaaddpath(path_to_fusion);
end 
%% Load data
mode_train = 1;
mode_test = 1;
RTC_delai = 1;
case_path = "../results/dcao/";
slopes_cl = fitsread(case_path+'integrator/saxoplus_KL_res.fits');
command_cl = fitsread(case_path+'integrator/saxoplus_KL_u.fits');

dist_matrix = command_cl(1:end-RTC_delai,:)+slopes_cl(1+RTC_delai:end,:);

%%
n_modes = 400;
dummy_vect = [1,2,3,4,5,6,50,100,200,300,n_modes+1];
n_controllers = length(dummy_vect);

bandwidth = [1,1,1,1,1,1,1,1,1,1]*300/15;
max_control_gain = [1,1,1,1,1,1,1,1,1,1];
% order = [3,3,3,3,3,3,3,3,3,3];
order = [1,1,1,1,1,1,1,1,1,1];
%%
fs = 300;


%%

fft_size = 200;
[psd_mat,f] = compute_psd_welch(dist_matrix,fft_size,fs);

% psd_mat(20:end,:) = repmat(psd_mat(20,:),size(psd_mat(20:end,:),1),1);

%%

% figure()
% semilogx(f,10*log10(psd_mat))
% psd(1:4) = psd(5);
% G = tf([1],[1,0,0],1/fs); % 2 samples delay
G = tf([1],[1,0],1/fs); % 1 samples delay

%% Load data
max_order =  max(order);


Kdd_matrix = zeros(max_order+1,2,n_modes);
Kdd_numerator = zeros(max_order+1,n_modes);
Kdd_denominator = zeros(max_order+1,n_modes);
tic
for i = 1:n_controllers-1
    psd_mat_max = max(psd_mat(:,dummy_vect(i):dummy_vect(i+1)-1),[],2);
    w = f*2*pi;
    W1 = frd(psd_mat_max,w,1/fs);
    val = freqresp(W1, bandwidth(i)*2*pi);
    W1 = W1/val;
    W3 = tf(max_control_gain(i),1,1/fs);
    Kdd =  ao_dd_controller(fs,w,order(i),W1,W3,[],'fusion');
    Kdd_numerator(1:order(i)+1,dummy_vect(i):dummy_vect(i+1)-1) = repmat(cell2mat(Kdd.Numerator)',[1,dummy_vect(i+1)-dummy_vect(i)]);
    Kdd_denominator(1:order(i)+1,dummy_vect(i):dummy_vect(i+1)-1) = repmat(cell2mat(Kdd.Denominator)',[1,dummy_vect(i+1)-dummy_vect(i)]);
end
toc
Kdd_matrix(:,1,:) = Kdd_numerator;
Kdd_matrix(:,2,:) = Kdd_denominator;

%% Simulation
mode_test = 300;
g = 0.4;

K0 = tf([g,0],[1,-1],1/fs);
K = tf(Kdd_matrix(:,1,mode_test)',Kdd_matrix(:,2,mode_test)',1/fs);
sys_dd = feedback(1,G*K);
sys_int = feedback(1,G*K0);

t = 0:1/fs:size(dist_matrix,1)/fs-1/fs;

[ydd,~] = lsim(sys_dd,dist_matrix(:,mode_test),t);
[yint,~] = lsim(sys_int,dist_matrix(:,mode_test),t);

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

S_dd = feedback(1,G*Kdd);
dummy = S_dd*W1*val;
figure()
bodemag(dummy,S_dd,W1)

S_int = feedback(1,G*K0);
dummy = S_int*W1*val;
figure()
bodemag(dummy,S_int,W1)

%% Save controller coefficients
% fitswrite(Kdd_matrix,'Kdd_matrix.fits')
Kdd_matrix = reshape(Kdd_matrix,(max_order+1)*2,n_modes);
% Kdd_matrix(max_order+2,:) = 0;
fitswrite(Kdd_matrix,'Kdd_matrix.fits')
% % Kdd = circshift(Kdd,-2,2); % put tip tilt controllers at the end 
% save('../Kdd_ProxCen','Kdd_matrix');