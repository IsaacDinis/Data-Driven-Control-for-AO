case_path = "results/with_defects/";
mode = 1;
% phase_ol = fitsread(case_path+'ol/saxoplus_KL_phase_res.fits');
% phase_cl = fitsread(case_path+'saxoplus_KL_phase_res.fits');
slopes_cl = fitsread(case_path+'HODM_res.fits');
command_cl = fitsread(case_path+'command.fits');

% phase_ol = phase_ol(:,mode);
% phase_cl = phase_cl(:,mode);
slopes_cl = slopes_cl(:,mode);
command_cl = command_cl(:,mode);
%%
fs = 4000;
gain = 0.8; % integrator gain
RTC_delai = 1; % number of RTC delay frames 
% delay = 1.83
K = tf([gain,0],[1,-0.99],1/fs); % integrator transfer function

dummy = zeros(1,1+RTC_delai);
dummy(1) = 1;
RTC = tf([1],dummy,1/fs); % RTC transfer function

WFS = tf([1,1],[2,0],1/fs); % WFS transfer function


% sys_ol = WFS*RTC*K; % open loop transfer function
sys_ol = RTC*K; % open loop transfer function
sys_cl = feedback(1,sys_ol); % closed loop transfer function

%%
sim_data = command_cl(1:end-RTC_delai)+slopes_cl(1+RTC_delai:end);
t = 0:1/fs:size(sim_data,1)/fs-1/fs;
sim_data = sim_data%+sin(2*pi*974*t)';



[res_sim,~] = lsim(sys_cl,sim_data,t);
res_sim = res_sim(20:end);
%%
fft_size = 500;

[psd_slopes_cl_welch,f] = compute_psd_welch(slopes_cl,fft_size,fs);
[psd_slopes_cl_fft,f] = compute_psd_fft(slopes_cl,fft_size,fs);

[psd_phase_ol_welch,f] = compute_psd_welch(sim_data,fft_size,fs);
[psd_phase_ol_fft,f] = compute_psd_fft(sim_data,fft_size,fs);

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
semilogx(f,20*log10(psd_phase_ol_fft))
hold on 
semilogx(f,20*log10(psd_res_sim_fft))
semilogx(f,-20*log10(psd_phase_ol_fft)+20*log10(psd_res_sim_fft))
title("PSD")
ylabel('mag [dB]')
xlabel('freq [Hz]')
legend('disturbance','residual','rejection')

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

%%
% fft_size = 2000;
% [psd_res_sim_fft,f] = compute_psd_fft(sim_data,fft_size,fs);
% [psd_slopes_cl_fft,f] = compute_psd_fft(phase_ol,fft_size,fs);
% figure()
% semilogx(f,20*log10(psd_res_sim_fft))
% hold on 
% semilogx(f,20*log10(psd_slopes_cl_fft))
% title("KL mode residual PSD, 0 frame delay in matlab, no WFS response")
% legend('matlab','compass')