import numpy as np
import astropy.io.fits as pfits
from matplotlib import pyplot as plt
from scipy.linalg import circulant

def compute_psd_welch(data, fft_size, fs):
    if data.ndim == 1:
        data = data[:, np.newaxis]

    n_modes = data.shape[1]

    if fft_size % 2 == 0:
        fft_size += 1

    window_size = fft_size * 2 - 1

    n_frames = (data.shape[0] // window_size) * 2 - 1
    spectrogram = np.zeros((fft_size, n_frames, n_modes))

    window = np.hamming(window_size) * 2

    for mode in range(n_modes):
        for i in range(n_frames):
            data_w = data[(i * fft_size):(i * fft_size + window_size), mode]
            data_w = data_w * window
            psd_w = np.abs(np.fft.fft(data_w)) / window_size
            psd_w = psd_w[:fft_size]
            psd_w[1:] = 2 * psd_w[1:]
            spectrogram[:, i, mode] = psd_w

    spectrogram = spectrogram.squeeze()
    psd = np.mean(spectrogram, axis=1).squeeze()

    f = fs * np.arange(0, (window_size // 2) + 1) / window_size

    return psd, f, spectrogram



fs = 1000
order = 3
train_size = 10000



t = np.arange(0,2*train_size/fs,1/fs)

# dist = np.sin(t*2*np.pi)

# dist = pfits.getdata('tilt_dist.fits').squeeze()
dist = pfits.getdata('P.fits').squeeze()
plt.figure()
plt.plot(dist)
PHI = np.zeros((train_size,order))


for i in range(order):
    PHI[:,i] = np.roll(dist,-i)[:train_size]
ar_coef = np.linalg.solve(PHI.T@PHI,PHI.T@dist[order:train_size+order])
ar_coef = ar_coef[::-1]


x_pred_1step = np.zeros_like(dist)
x_pred_2step = np.zeros_like(dist)

for i in range(order, dist.shape[0]):
    x_pred_1step[i] = np.dot(ar_coef, dist[i:i - order:-1])
    x_pred_2step[i] = np.dot(ar_coef[1:], dist[i:i - order + 1:-1])+x_pred_1step[i]*ar_coef[0]

end_plot = order+1000
plt.figure()
plt.plot(x_pred_2step[order:end_plot])
plt.plot(dist[order:end_plot])
plt.plot(dist[order+2:end_plot+2])
plt.legend(('pred','meas','true'))
plt.show()

print(np.std(dist[order+2:]-x_pred_2step[order:-2]))
print(np.std(dist[order+2:]-dist[order:-2]))

# psd_res,f,_ = compute_psd_welch(dist[order+2:]-x_pred_2step[order:-2], 300, fs)
# psd_res,f,_ = compute_psd_welch(dist, 300, fs)
# plt.figure()
# plt.semilogx(f, 20 * np.log10(psd_res))
# plt.title('residual PSD')
# plt.xlabel("frequency [Hz]")
# plt.ylabel("magnitude [dB]")
# plt.grid()
# plt.show()

plt.figure()
plt.plot(ar_coef)
plt.show()