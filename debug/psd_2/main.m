case_path = "results/bright1/04_9/no_mod/";
mode = 1;
phase_ol = fitsread(case_path+'ol/saxoplus_KL_phase_res.fits');
phase_cl = fitsread(case_path+'cl_g02_2frames/saxoplus_KL_phase_res.fits');
slopes_cl = fitsread(case_path+'cl_g02_2frames/saxoplus_KL_res.fits');
command_cl = fitsread(case_path+'cl_g02_2frames/saxoplus_KL_u.fits');

phase_ol = phase_ol(:,mode);
phase_cl = phase_cl(:,mode);
slopes_cl = slopes_cl(:,mode);
command_cl = command_cl(:,mode);
%%
fs = 2760;
gain = 0.6; % integrator gain
RTC_delai = 2; % number of RTC delay frames 

K = tf([gain,0],[1,-1],1/fs); % integrator transfer function

dummy = zeros(1,1+RTC_delai);
dummy(1) = 1;
RTC = tf([1],dummy,1/fs); % RTC transfer function

WFS = tf([1,1],[2,0],1/fs); % WFS transfer function


% sys_ol = WFS*RTC*K; % open loop transfer function
sys_ol = RTC*K; % open loop transfer function
sys_cl = feedback(1,sys_ol); % closed loop transfer function
%%
t = 0:1/fs:size(command_cl,1)/fs-1/fs;
% [res_sim,~] = lsim(sys_cl,command_cl+slopes_cl,t);
[res_sim,~] = lsim(sys_cl,phase_ol,t);
res_sim = res_sim(20:end);
%%
fft_size = 300;

[psd_slopes_cl_welch,f] = compute_psd_welch(slopes_cl,fft_size,fs);
[psd_slopes_cl_fft,f] = compute_psd_fft(slopes_cl,fft_size,fs);

[psd_phase_ol_welch,f] = compute_psd_welch(phase_ol,fft_size,fs);
[psd_phase_ol_fft,f] = compute_psd_fft(phase_ol,fft_size,fs);

[psd_res_sim_welch,f] = compute_psd_welch(res_sim,fft_size,fs);
[psd_res_sim_fft,f] = compute_psd_fft(res_sim,fft_size,fs);
%%
% figure()
% semilogx(f,20*log10(psd_res_sim_welch))
% hold on 
% semilogx(f,20*log10(psd_slopes_cl_welch))
% title("KL mode residual PSD")
% legend('matlab','compass')


figure()
semilogx(f,20*log10(psd_res_sim_fft))
hold on 
semilogx(f,20*log10(psd_slopes_cl_fft))
title("KL mode residual PSD, 0 frame delay in matlab, no WFS response")
legend('matlab','compass')

figure()
plot(res_sim)
hold on;
plot(slopes_cl)

