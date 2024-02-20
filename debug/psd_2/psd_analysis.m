case_path = "results/bright1/04_9/no_mod/";
phase_ol = fitsread(case_path+'ol/saxoplus_KL_phase_res.fits');
slopes_ol = fitsread(case_path+'ol/saxoplus_KL_res.fits');
phase_cl = fitsread(case_path+'cl_no_noise/saxoplus_KL_phase_res.fits');
slopes_cl = fitsread(case_path+'cl_no_noise/saxoplus_KL_res.fits');
command_cl = fitsread(case_path+'cl_no_noise/saxoplus_KL_u.fits');

mode = 10;

figure()
plot(phase_ol(:,mode))
hold on;
plot(slopes_ol(:,mode))
plot(command_cl(:,mode))
legend('phase','slopes','command')

figure()
plot(slopes_cl(:,mode))
hold on;
plot(phase_cl(:,mode))
plot(phase_ol(:,mode)-command_cl(:,mode))

legend('slopes','phase cl','phase ol')

figure()
plot(phase_cl(:,mode)-phase_ol(:,mode)+command_cl(:,mode))

%% sanity check 
fs = 1000;
T = 5;
t = 0:1/fs:T-1/fs;
fft_size = 1000;
s = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t)+4;
s = s';
[psd,f] = compute_psd_welch(s,fft_size,fs);
figure()
plot(f,psd)
sum(psd)
%%
fs = 2760;
fft_size = 1000;
mode = 1;

[psd_slopes_cl_welch,f] = compute_psd_welch(slopes_cl,fft_size,fs);
[psd_slopes_cl_fft,f] = compute_psd_fft(slopes_cl,fft_size,fs);
[psd_phase_ol_welch,f] = compute_psd_welch(phase_ol,fft_size,fs);
[psd_phase_ol_fft,f] = compute_psd_fft(phase_ol,fft_size,fs);

figure()
semilogx(f,20*log10(psd_slopes_cl_welch(:,mode)))
hold on 
semilogx(f,20*log10(psd_slopes_cl_fft(:,mode)))

figure()
semilogx(f,20*log10(psd_phase_ol_welch(:,mode)))
hold on 
semilogx(f,20*log10(psd_phase_ol_fft(:,mode)))


figure()
semilogx(f,20*log10(psd_slopes_cl_welch(:,mode))-20*log10(psd_phase_ol_welch(:,mode)))
hold on 
semilogx(f,20*log10(psd_slopes_cl_fft(:,mode))-20*log10(psd_phase_ol_fft(:,mode)))



