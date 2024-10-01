import numpy as np
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