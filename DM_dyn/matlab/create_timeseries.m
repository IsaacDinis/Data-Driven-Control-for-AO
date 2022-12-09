dist_matrix_path = '../data/single_mode_dist_ProxCen_1_4000.mat';
dist = load(dist_matrix_path).data';
dist = double(dist(:));
fs = 2000;
t = 0:1/fs:size(dist,1)/fs-1/fs;
plop = timeseries(dist, t);
save('../simulink/single_mode_dist_ProxCen_1_4000_ts.mat','plop','-v7.3')

% save 'dist.mat' -v7.3 timeseries_object
