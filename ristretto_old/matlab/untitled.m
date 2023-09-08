dist_matrix_path = '../data/single_mode_dist_ProxCen_1.mat';
dist_matrix = load(dist_matrix_path).data';
dist_matrix = dist_matrix(2:end);

fs = 1000;
bandwidth = 400;

window_size = 49*5;
n_average = 5*2*20;
hann_window = hann(length(dist_matrix));

[psd, f] = compute_psd(dist_matrix,n_average,window_size,fs);
w = f*2*pi;

% [psd,w] = periodogram(dist_matrix,hann_window,1000,fs);

W1 = frd(psd,w,1/fs);
val = freqresp(W1, bandwidth*2*pi);
% val = interp1(f,psd,bandwidth);
W1 = W1/val;

W12 = makeweight(10^(54.5/20),[152,10^(34/20)],10^(20/20),1/fs,2);

figure()
bodemag(W1,W12)


% pmtm(dist_matrix,{7.5,131},(w/fs),fs)