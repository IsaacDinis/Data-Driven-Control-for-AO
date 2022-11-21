%%
dist_ProxCen = load('../data/single_mode_dist_ProxCen_1_2000.mat').data';
dist_ProxCen = dist_ProxCen(2:end);

dist_PDS70 = load('../data/single_mode_dist_PDS70_1_2000.mat').data';
dist_PDS70 = dist_PDS70(2:end);

dist_nonoise = load('../data/single_mode_dist_nonoise_1_2000.mat').data';
dist_nonoise = dist_nonoise(2:end);
%%
window_size = 999;
n_average = 20;
fs = 2000;

[psd_ProxCen, f] = compute_psd(dist_ProxCen,n_average,window_size,fs);
[psd_PDS70, f] = compute_psd(dist_PDS70,n_average,window_size,fs);
[psd_nonoise, f] = compute_psd(dist_nonoise,n_average,window_size,fs);

psd_ProxCen = psd_ProxCen/max(psd_ProxCen);
psd_PDS70 = psd_PDS70/max(psd_PDS70);
psd_nonoise = psd_nonoise/max(psd_nonoise);
%%

figure()
semilogx(f,10*log10(psd_nonoise))
hold on;
semilogx(f,10*log10(psd_ProxCen))
semilogx(f,10*log10(psd_PDS70))
title('Disurbance PSD')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('no photon noise','Prox Cen.','PDS70','Interpreter','latex')
grid()
make_it_nicer()

set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
export_fig ../plot/psd.pdf -transparent