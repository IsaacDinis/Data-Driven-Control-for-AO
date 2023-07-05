import numpy as np
from scipy.signal import correlate
from matplotlib import pyplot as plt


def autocorr(u):
    N = u.size
    R = np.zeros(2 * N - 1)
    y = np.concatenate((u,u,u))
    for h in range(-N+1,N - 1):
        sum = 0
        for k in range(N):
            sum += u[k] * y[k - h + N]
        R[h + N] = sum
    h = np.arange(-N,N - 1)
    R = R / N
    return R,h


def compute_psd_fft(data, n_average, window_size,fs):
    try:
        n_modes = data.shape[1]
    except IndexError:
        n_modes = 1
        data = data.reshape((data.shape[0],1))
    hann_window = np.hanning(window_size )
    psd = np.zeros((window_size, n_modes))

    for mode in range(n_modes):
        for i in range(n_average):
            data_w = data[i*window_size:(i+1)*window_size, mode] * hann_window # TODO check if window useful
            psd[:, mode] = psd[:, mode] + np.abs((np.fft.fft(data_w)/data_w.size))**2

    f = fs / (2 * window_size + 1) * np.arange(0,2 * window_size)

    # remove frequencies above nyquist
    f = f[:window_size]
    psd = psd[:window_size,:]

    # freq = np.fft.fftfreq(psd.size,1/fs)
    # psd = psd[:int(psd.size / 2)]
    # freq = freq[:int(freq.size/2)]

    return psd, freq


def compute_psd_xcorr(data, n_average, window_size, fs):
    try:
        n_modes = data.shape[1]
    except IndexError:
        n_modes = 1
        data = data.reshape((data.shape[0],1))
    hann_window = np.hanning(2 * window_size+1)
    psd = np.zeros((2 * window_size + 1, n_modes))

    for mode in range(n_modes):
        for i in range(n_average):
            data_xcorr = correlate(data[i*window_size:(i+1)*window_size+1, mode],data[i*window_size:(i+1)*window_size+1, mode])/window_size
            # data_xcorr,_ = autocorr(data[i * window_size:(i + 1) * window_size + 1, mode]) # more accurate but slower
            data_xcorr_w = hann_window * data_xcorr
            psd[:, mode] = psd[:, mode] + np.abs(np.fft.fft(data_xcorr_w)/window_size)

    psd = psd / n_average
    f = fs / (2 * window_size + 1) * np.arange(0,2 * window_size)

    # remove frequencies above nyquist
    # f = f[:window_size]
    # psd = psd[:window_size,:]

    return psd, f


if __name__ == '__main__':
    T = 50
    fs = 1000
    t = np.arange(0,T,1/fs)
    signal = 3*np.sin(5*2*np.pi*t) + 7*np.sin(100*2*np.pi*t) + 2*np.sin(367*2*np.pi*t)

    psd, freq = compute_psd_xcorr(signal, 1, 50 * fs - 1, fs)
    psd2, freq2 = compute_psd_fft(signal, 1, 50 * fs - 1, fs)

    # corr,h = autocorr(a)

    # plt.plot(h, corr)
    # plt.show()

    print(np.std(signal) ** 2)
    print(np.sum(psd))
    print(np.sum(psd2))

    plt.plot(freq, psd)
    plt.show()
    plt.plot(freq2, psd2)
    plt.show()