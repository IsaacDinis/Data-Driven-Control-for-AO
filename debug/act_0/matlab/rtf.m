ol = fitsread('../compass/results/ol/zernike_res.fits');
cl = fitsread('../compass/results/cl/zernike_res.fits');
n_average = 10;
window_size = 390;
fs = 2000;

[psd_average_ol,f,psd_std_ol]=  compute_psd(ol,n_average,window_size,fs);
psd_average_ol = psd_average_ol(2:end,:);
f = f(2:end);
psd_std_ol = psd_std_ol(2:end,:);

[psd_average_cl,f,psd_std_cl]=  compute_psd(cl,n_average,window_size,fs);
psd_average_cl = psd_average_cl(2:end,:);
f = f(2:end);
psd_std_cl = psd_std_cl(2:end,:);

[psd_average_rtf,f_rtf,psd_std_rtf]=  compute_psd(cl(2:end,1)./ol(2:end,1),n_average,window_size,fs);
psd_average_rtf = psd_average_rtf(2:end,:);

psd_std_cl = psd_std_cl(2:end,:);

figure()
plot(ol(:,1))


figure()
semilogx(f,10*log10(psd_average_ol(:,1)))
hold on;
semilogx(f,10*log10(psd_average_cl(:,1)))
% semilogx(f,10*log10(psd_average_cl(:,1))-10*log10(psd_average_ol(:,1)))
semilogx(f,10*log10(psd_average_rtf))
title('Tilt PSD')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('open-loop','closed-loop','RTF','Interpreter','latex')
make_it_nicer()


figure()
semilogx(f,10*log10(psd_average_cl(:,1))-10*log10(psd_average_ol(:,1)))

title('First KL dcao')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()

%%

K = tf([0.5],[1,-1],1/fs);
G = tf([1],[1,0,0],1/fs);
S = feedback(1,G*K);
K2 = tf([0.5],[1,-0.99],1/fs);
S2 = feedback(1,G*K2);
figure()
bodemag(S,S2)
title('RTF')
legend('integrateur standard','leak 1%')

[mag,phase,wout] = bode(S,f*2*pi);
mag = squeeze(mag);
% figure()
% bodemag(K,K3)

figure()
semilogx(wout,10*log10(mag))

figure()
semilogx(f,10*log10(psd_average_ol(:,1))+10*log10(mag))

title('First KL dcao')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
make_it_nicer()