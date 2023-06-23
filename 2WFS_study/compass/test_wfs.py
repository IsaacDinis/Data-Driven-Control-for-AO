from hcipy.field import make_pupil_grid
from hcipy.mode_basis import make_zernike_basis
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

def compute_psf(phase, mask, step=1, pad=0):

    # phase[mask > 0] -= np.mean(phase[mask > 0])

    e_field = np.exp(phase*1j)*mask
    e_field /= np.sum(mask)
    e_field = e_field[::step, ::step]
    e_field = np.pad(e_field, pad, 'constant', constant_values=0)

    # plt.imshow(np.abs(e_field))
    # plt.colorbar()
    # plt.show()

    psf = np.fft.fftshift(np.abs(np.fft.fft2(e_field))**2)

    return psf

def plot_psf(psf,pupil_diam, step=1, pad=0):

    fs = pupil_diam/step
    res = fs/((pupil_diam/step)+2*pad)
    ticks = np.arange(-fs/2, fs/2, res)

    plt.plot(ticks,psf[:, int(psf.shape[0] / 2)])
    # plt.yscale('log')
    plt.show()

    plt.imshow(psf,extent=[-fs,fs,-fs,fs])
    plt.colorbar()
    plt.show()

if __name__ == '__main__':
    print('hello')
    # lambd = 1.65

    # D = 8
    pupil_diam = 128
    pad = 0
    step = 1

    pupil_grid = make_pupil_grid(pupil_diam)
    zernike_basis = make_zernike_basis(3, 1, pupil_grid)
    tilt = zernike_basis[2].shaped
    tilt = np.array(tilt)
    pupil_valid = np.array(zernike_basis[0].shaped)
    # pupil_valid = np.array(pupil_valid)*1
    # plt.imshow(pupil_valid)
    # plt.colorbar()
    # plt.show()
    tilt *= np.pi/np.max(tilt)
    psf = compute_psf(tilt, pupil_valid, step=step, pad=pad)
    plot_psf(psf, pupil_diam, step=step, pad=pad)

    fs = pupil_diam/step
    res = fs/((pupil_diam/step)+2*pad)
    ticks = np.arange(-fs, fs, res*2)
