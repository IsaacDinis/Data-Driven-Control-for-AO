P2M_DM0 = fitsread('../compass/calib_mat/P2M_DM0.fits');
P2M_DM1 = fitsread('../compass/calib_mat/P2M_DM1.fits');

M2P_DM0 = fitsread('../compass/calib_mat/M2P_DM0.fits');
M2P_DM1 = fitsread('../compass/calib_mat/M2P_DM1.fits');

M2V_DM0 = fitsread('../compass/calib_mat/M2V_DM0.fits');
phase = fitsread('../compass/calib_mat/target_phase.fits');
tip = fitsread('../compass/calib_mat/target_phase_tip.fits');
% 
% phase = ones(11304,1);
% 
% M_DM0 = P2M_DM0*phase;
% M_DM1 = P2M_DM1*phase;
% 
% [V,D] = eig(M2P_DM0'*M2P_DM0);
% 
% figure()
% imagesc(U)

% M2V_DM0 = M2V_DM0./max(M2V_DM0,[],'all');
% volt = ones(88,1);
% sum(M2V_DM0(:,3).*volt)

% figure()
% plot(M2V_DM0(:,1))
% 
% figure()
% plot(M2P_DM0(:,2))
% 
% figure()
% plot(tip)
tip =tip';
dm = M2P_DM1(:,1)./max(M2P_DM1(:,1)).*max(tip);
dm = dm*1;
figure()
plot(tip)
figure()
plot(dm)
result = tip+dm;
std(result)
std(tip)
figure()
plot(result)