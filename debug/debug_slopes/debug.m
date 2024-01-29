
slopes = fitsread('dcao_2/saxoplus_KL_res.fits');
phase = fitsread('dcao_2/saxoplus_KL_phase_res.fits');
phase_std = fitsread('dcao_2/saxoplus_std_phase.fits');
phase_std = phase_std.*phase_std;
phase_rms = rms(phase,2)/sqrt(121140);
slopes_rms = rms(slopes,2)/sqrt(8064);
mode = 1;
figure()
plot(slopes_rms)
hold on;
plot(phase_rms)
plot(phase_std)
legend('slope','kl phase','phase')

mean(phase_rms./phase_std)
mean(phase_std)


% figure()
% plot(slopes(:,mode)/(8064))
% hold on;
% plot(phase(:,mode)/(121140))

figure()
plot(slopes(:,mode))
hold on;
plot(phase(:,mode))

%%
atm = fitsread('dcao_2/atm_KL_phase.fits');
command = fitsread('dcao_2/saxoplus_KL_u.fits');
slopes_ol = fitsread('dcao_2/saxoplus_KL_ol_res.fits');
figure()
plot(atm(:,mode))
hold on;
plot(command(:,mode))
plot(slopes_ol(:,mode))

figure()
plot(slopes(:,mode))
hold on;
plot(phase(:,mode))

figure()
plot(atm(1:end-2,mode))
hold on;
% plot(command(3:end,mode))
% plot(command(3:end,mode)+slopes(1:end-2,mode))
plot(command(3:end,mode)+phase(1:end-2,mode))
% plot(slopes_ol(1:end-2,mode))
%%
fs = 1000;
G = tf([1],[1,0,0],1/fs); % 1 samples delay
K0 = tf([0.5,0],[1,-1],1/fs);
S = feedback(1,G*K0);
figure()
step(S)

%%
mode = 1;


fs = 2760;



window_size = 200;
n_average = 60;

[psd_atm, f] = compute_psd(atm(:,mode),n_average,window_size,fs);
[psd_slopes, f] = compute_psd(slopes_ol(:,mode),n_average,window_size,fs);
[psd_command, f] = compute_psd(command(:,mode),n_average,window_size,fs);

figure()
semilogx(10*log10(psd_atm))
hold on;
semilogx(10*log10(psd_slopes))
semilogx(10*log10(psd_command))
legend('atm','slopes','command')
% psd(25:end) = psd(25);
