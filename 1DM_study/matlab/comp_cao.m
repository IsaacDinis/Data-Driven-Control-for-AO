command_CAO = fitsread('../data/command_CAO.fits');
command_dCAO = fitsread('../data/command_dCAO_volt.fits');
res_dCAO = fitsread('../data/res_dCAO_volt.fits');
res_CAO = fitsread('../data/res_CAO.fits');
fs = 4000;
t = 0:1/fs:0.05-1/fs;
t_start = 1001;
t_end = 1200;

%%
figure()

plot(t,command_CAO(t_start:t_end))
hold on
plot(t,command_dCAO(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Command 1st KL')
legend('CAO','dCAO','Interpreter','latex', 'Location','southwest');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
%%
figure()

plot(t,res_CAO(t_start:t_end))
hold on
plot(t,res_dCAO(t_start:t_end))
xlabel('Time (s)')
ylabel('Res. amp.')
title('Rsidual 1st KL')
legend('CAO','dCAO','Interpreter','latex', 'Location','southwest');
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
make_it_nicer()
%%
n = 2;
w = 100;

[psd,f] = compute_psd((res_dCAO)',n,w,fs);

[inter_psd,f] = compute_psd((res_CAO)',n,w,fs);

figure()
plot(f(1:end),psd(1:end))
hold on;
plot(f(1:end),inter_psd(1:end))


%% PSD
n = 13;
w = 300;

[psd_CAO,f] = compute_psd(res_CAO',n,w,fs);

[psd_dCAO,f] = compute_psd(res_dCAO',n,w,fs);

sum_psd_CAO = sum(psd_CAO);
sum_psd_dCAO = sum(psd_dCAO);
figure()
plot(f(1:167),psd_CAO(1:167))
hold on;
plot(f(1:167),psd_dCAO(1:167))
legend('CAO','dCAO','Interpreter','latex');
title('residual 1st Kl PSD, using getcommand')
xlabel('Frequency (Hz)')
ylabel('Amp.')
% xlim([0,2050])
make_it_nicer()
set(gcf, 'Position',  [100, 100, 700, 450])
set(gcf,'PaperType','A4')
% export_fig ../plot/case_1_psd.pdf -transparent
