mode = 1;
dist_matrix_path = '../results/bright1_03_8ms/dcao_ol/saxoplus_KL_res.fits';
% dist_matrix_path = '../data/single_mode_dist.mat';
dist_matrix = fitsread(dist_matrix_path);
dist_matrix = dist_matrix(:,mode);
fs = 2760;
figure()
[s,f,t] = stft(dist_matrix,fs,Window=kaiser(256,5),OverlapLength=220,FFTLength=512);
sdb = mag2db(abs(s));
mesh(t,f/1000,sdb);

cc = max(sdb(:))+[-60 0];
ax = gca;
set(gca, 'YScale', 'log')
ax.CLim = cc;
view(2)
colorbar


figure()
periodogram(dist_matrix)
set(gca, 'XScale', 'log')