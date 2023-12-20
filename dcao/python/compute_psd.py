import numpy as np
from scipy.signal import correlate
from matplotlib import pyplot as plt



def compute_psd_fft(data, n_average, window_size,fs):
    try:
        n_modes = data.shape[1]
    except IndexError:
        n_modes = 1
        data = data.reshape((data.shape[0],1))

    psd = np.zeros((window_size, n_modes))

    for mode in range(n_modes):
        for i in range(n_average):
            data_w = data[i*window_size:(i+1)*window_size, mode]
            psd[:, mode] += np.abs((np.fft.fft(data_w)))**2/data_w.size

    freq = np.fft.fftfreq(psd.shape[0],1/fs)

    # remove frequencies above nyquist
    psd = psd[:int(psd.shape[0] / 2),:]
    freq = freq[:int(freq.size/2)]

    psd *= 2/(n_average*fs) # normalize energy because of removing negative parts and number of averages
    return psd, freq

def compute_psd_welch(data, N, R,fs):
    try:
        n_modes = data.shape[1]
    except IndexError:
        n_modes = 1
        data = data.reshape((data.shape[0],1))

    psd = np.zeros((window_size, n_modes))

    for mode in range(n_modes):
        for i in range(n_average):
            data_w = data[i*window_size:(i+1)*window_size, mode]
            psd[:, mode] += np.abs((np.fft.fft(data_w)))**2/data_w.size

    freq = np.fft.fftfreq(psd.shape[0],1/fs)

    # remove frequencies above nyquist
    psd = psd[:int(psd.shape[0] / 2),:]
    freq = freq[:int(freq.size/2)]

    psd *= 2/(n_average*fs) # normalize energy because of removing negative parts and number of averages
    return psd, freq



if __name__ == '__main__':
    n_average = 10 # to decrase variance
    window_size = 1000 # defines number of frequency bins between 0 and nyquist
    # /!\ window_size * n_average < length of signal

    T = 50
    fs = 1000
    t = np.arange(0,T,1/fs)
    signal = 3*np.sin(5*2*np.pi*t) + 7*np.sin(100*2*np.pi*t) + 2*np.sin(367*2*np.pi*t)


    psd_fft, freq_fft = compute_psd_fft(signal, n_average, window_size, fs)

    print(np.std(signal) ** 2)

    print(np.sum(psd_fft)*freq_fft[1])

    plt.figure()
    plt.plot(freq_fft, psd_fft)
    plt.show()