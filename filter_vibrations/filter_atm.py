import numpy as np
import astropy.io.fits as pfits
from matplotlib import pyplot as plt
import dd4ao
from dd_utils import *

traj = pfits.getdata('traj/traj_pol_0_tt.fits').T
traj = traj[:,0]
fs = 1380
n_fft = 10000
controller_order = 10
fit_order = 5
n_times = 50

# plt.figure()
# plt.plot(traj)
# plt.title('pseudo open loop tilt')
# plt.show()

traj_psd,f, _  = compute_fft_mag_welch(traj,n_fft,fs)
f = f[1:]
traj_psd = traj_psd[1:]
# coefficients = np.polyfit(np.arange(f.shape[0]), np.log10(traj_psd), fit_order)
# polynomial = np.poly1d(coefficients)
# y_fit = polynomial(np.arange(f.shape[0]))

# plt.figure()
# plt.semilogx(np.arange(f.shape[0]), np.log10(traj_psd))
# plt.semilogx(np.arange(f.shape[0]),y_fit)

w_log,traj_psd_log = interp_log(f*2*np.pi,traj_psd,300)

mask = ~((w_log >= 26) & (w_log <= 720))

w_log_mask = np.concatenate((w_log[mask],np.repeat(w_log[68], n_times*100), np.repeat(w_log[152], n_times), np.repeat(w_log[194], n_times)))
traj_psd_log_mask = np.concatenate((traj_psd_log[mask], np.repeat(traj_psd_log[68], n_times*100), np.repeat(traj_psd_log[152], n_times), np.repeat(traj_psd_log[194], n_times)))


coefficients = np.polyfit(np.log10(w_log_mask), np.log10(traj_psd_log_mask), fit_order)
polynomial = np.poly1d(coefficients)
y_fit = polynomial(np.log10(w_log))
y_fit_2 = y_fit.copy()
y_fit_2[:68] = np.log10(traj_psd_log[:68])
# y_fit_2[245:] = -5
# y_fit_2[245:] = -1000
y_fit_2[245:] = -5
plt.figure()
plt.semilogx(w_log, y_fit)
plt.semilogx(w_log, np.log10(traj_psd_log))

plt.figure()
plt.semilogx(w_log, 10**y_fit)
plt.semilogx(w_log, traj_psd_log)
plt.semilogx(w_log, 10**y_fit_2)

# w_log,traj_psd_log = interp_log(f*2*np.pi,traj_psd,300)

# # mask = ~((w_log >= 4*2*np.pi) & (w_log <= 1122*np.pi))

# # # Apply the mask to both arrays
# # w_log = w_log[mask]
# # traj_psd_log = traj_psd_log[mask]

G = G_tf(1, fs)
G_resp = freqresp(G, w_log)

K0_num = np.array([0.2,0])
K0_den = np.array([1,-1])
K0 = ct.tf(K0_num,K0_den, 1/fs)


bandwidth = 17

dd = dd4ao.DD4AO(w_log, G_resp, 1/y_fit_2, controller_order, bandwidth, fs, K0_num, K0_den)
dd.compute_controller()

K = dd.K

res_K0, u_K0 = evaluate_K_performance(K0,G,traj,fs)
res_K, u_K = evaluate_K_performance(K,G,traj,fs)
# plot_combined(G, K, K0, traj_psd, f, traj, res_K, 300, fs,0, bandwidth)
# plot_combined(G, K, K0, 10**y_fit_2, w_log, res_K0, res_K, 300, fs,0, bandwidth)


dist_psd = 10**y_fit_2

w_log = w_log[:-1]
dist_psd = dist_psd[:-1]
f = f[:-1]
traj_psd = traj_psd[:-1]

val = np.interp(115 * 2 * np.pi, w_log.squeeze(), dist_psd.squeeze())
dist_psd = dist_psd / val

K_cl = ct.feedback(1, G*K)
K_cl_freqresp = np.abs(freqresp(K_cl, f*2*np.pi))
plt.figure()
plt.semilogx(f, 20 * np.log10(K_cl_freqresp))
plt.semilogx(w_log/2/np.pi, 20 * np.log10(1/dist_psd))
plt.legend(('datadriven', 'disturbance^-1'))
# plt.set_title('sensitivity function')
# plt.set_xlabel("frequency [Hz]")
# plt.set_ylabel("magnitude [dB]")
plt.grid()


val = np.interp(127, f.squeeze(), traj_psd.squeeze())
traj_psd_2 = traj_psd / val
plt.figure()
plt.semilogx(f, 20 * np.log10(K_cl_freqresp))
plt.semilogx(f, 20 * np.log10(1/traj_psd_2))
plt.legend(('datadriven', 'disturbance^-1'))
# plt.set_title('sensitivity function')
# plt.set_xlabel("frequency [Hz]")
# plt.set_ylabel("magnitude [dB]")
plt.grid()

res_K_fft, _,_  = compute_fft_mag_welch(res_K,n_fft,fs)
res_K_fft = res_K_fft[1:-1]

plt.figure()
plt.semilogx(f, 20 * np.log10(traj_psd))
plt.semilogx(f, 20 * np.log10(res_K_fft))
plt.legend(('datadriven', 'disturbance^-1'))
# plt.set_title('sensitivity function')
# plt.set_xlabel("frequency [Hz]")
# plt.set_ylabel("magnitude [dB]")
plt.grid()


plt.figure()
plt.semilogx(f, (traj_psd))
plt.semilogx(f, (res_K_fft))
plt.legend(('datadriven', 'disturbance^-1'))
# plt.set_title('sensitivity function')
# plt.set_xlabel("frequency [Hz]")
# plt.set_ylabel("magnitude [dB]")
plt.grid()

# plt.figure()
# plt.semilogx(f, np.cumsum(traj_psd))
# plt.semilogx(f, np.cumsum(res_K_fft))
# plt.legend(('traj', 'vib'))
# # plt.set_title('sensitivity function')
# # plt.set_xlabel("frequency [Hz]")
# # plt.set_ylabel("magnitude [dB]")
# plt.grid()

