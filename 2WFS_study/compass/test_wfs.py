from hcipy.field import make_pupil_grid
from hcipy.mode_basis import make_zernike_basis
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

def compute_psf(phase, mask, step=1, pad=0):
    # phase *= phase_factor
    phase -= np.mean(phase[mask > 0])
    # print(np.exp(-np.std(phase[mask>0])**2))
    e_field = np.exp(phase*1j)*mask


    e_field = e_field[::step, ::step]
    e_field = np.pad(e_field, pad, 'constant', constant_values=0)

    # plt.imshow(np.angle(e_field))
    # plt.colorbar()
    # plt.show()


    psf = np.fft.fftshift(np.abs(np.fft.fft2(e_field))**2)

    plt.imshow(psf)
    plt.colorbar()
    plt.show()

    return psf

if __name__ == '__main__':
    print('hello')
    lambd = 1.65
    pupil_diam = 128
    pad = 64
    D = 8
    step = 1

    pupil_grid = make_pupil_grid(pupil_diam)
    zernike_basis = make_zernike_basis(3, 1, pupil_grid)
    tilt = zernike_basis[1].shaped
    pupil_valid = zernike_basis[0].shaped
    plt.imshow(pupil_valid)
    plt.colorbar()
    plt.show()
    tilt *= np.pi/np.max(tilt)
    psf = compute_psf(tilt, pupil_valid, step=step, pad=pad)

    fs = pupil_diam/step
    res = fs/((pupil_diam/step)+2*pad)
    ticks = np.arange(-fs / 2, fs / 2, res)
    plt.plot(ticks,psf[:, int(psf.shape[0] / 2)])
    plt.show()
    a = 0