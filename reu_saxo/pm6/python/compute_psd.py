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

def compute_psd_welch(data, window_size,fs):

    if window_size % 2 != 0:
        window_size += 1
    try:
        n_modes = data.shape[1]
    except IndexError:
        n_modes = 1
        data = data.reshape((data.shape[0],1))

    n_frames = int(np.floor(data.shape[0]/window_size)*2)-1
    psd = np.zeros((window_size, n_frames, n_modes))
    psd_average = np.zeros((window_size, n_modes))

    window = np.hamming(window_size)
    for mode in range(n_modes):
        for i in range(n_frames):
            data_w = data[i*int(window_size/2):i*int(window_size/2)+window_size, mode].copy()
            data_w *= window
            psd[:, i, mode] = np.abs((np.fft.fft(data_w)))**2/data_w.size

    freq = np.fft.fftfreq(psd.shape[0],1/fs)

    # remove frequencies above nyquist
    psd = psd[:int(psd.shape[0] / 2),:,:]
    freq = freq[:int(freq.size/2)]
    psd *= 2/fs # normalize energy because of removing negative parts and number of averages

    psd_average = np.mean(psd, axis=1)
    return psd_average, freq, psd



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

    fig0, ax0 = plt.subplots()
    ax0.plot(freq_fft, psd_fft)


    psd_fft, freq_fft, psd_map = compute_psd_welch(signal, window_size, fs)

    print(np.std(signal) ** 2)

    print(np.sum(psd_fft)*freq_fft[1])

    fig1, ax1 = plt.subplots()
    ax1.plot(freq_fft, psd_fft)


    fig2, ax2 = plt.subplots()
    ax2.imshow(psd_map,extent=[0,100,500,0])

    # ax2.set_yticks(list(range(len(freq_fft))))
    # ax2.set_yticks(range(len((freq_fft))),freq_fft)
    # ax2.set_yticks(range(len((freq_fft))))
    # ax2.set_yticklabels(freq_fft)

    plt.show()
