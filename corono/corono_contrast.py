import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

def compute_strehl(phase_turb,psf_diff_limit,mask):
    psf_turb = compute_psf(phase_turb,mask)
    strehl = np.max(psf_turb)/np.max(psf_diff_limit)
    return strehl

def compute_pupil_mask(pupil_diam, c_obs):
    center = (int(pupil_diam / 2), int(pupil_diam / 2))
    r1 = center[0]
    r2 = pupil_diam * c_obs / 2
    Y, X = np.ogrid[:pupil_diam, :pupil_diam]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = 1.*((dist_from_center <= r1) * (dist_from_center >= r2))
    return mask

def compute_psf(phase, mask, step=1, pad=0):
    phase -= np.mean(phase[mask > 0])

    e_field = np.exp(phase*1j)*mask
    e_field = e_field[::step, ::step]
    e_field = np.pad(e_field, pad, 'constant', constant_values=0)

    psf = np.fft.fftshift(np.abs(np.fft.fft2(e_field))**2)
    return psf


def apply_corono(phase,psf_diff_limit,mask,step=1,pad=0):
    phase -= np.mean(phase[mask > 0])

    e_field = np.exp(phase*1j)*mask

    strehl = compute_strehl(phase,psf_diff_limit,mask)

    e_field -= mask*np.sqrt(strehl)
    e_field = e_field[::step, ::step]
    e_field = np.pad(e_field, pad, 'constant', constant_values=0)

    psf_corono = np.fft.fftshift(np.abs(np.fft.fft2(e_field))**2)
    return psf_corono

def load_phase(path_to_phase,lamb):
    return np.genfromtxt(path_to_phase, delimiter=",")*2*np.pi/lamb

def average_corono(path_to_phase,path_to_phase_diff_limit,lamb,n_frames,n_pixels,c_obs,step = 1, pad = 0):
    phase_diff_limit = load_phase(path_to_phase_diff_limit,lamb)
    mask = compute_pupil_mask(n_pixels,c_obs)
    psf_diff_limit = compute_psf(phase_diff_limit, mask)
    psf_shape = int(n_pixels/step+2*pad)
    psf_corono = np.zeros((psf_shape,psf_shape))
    psf = np.zeros_like(psf_corono)

    for i in range(n_frames):
        phase_turb = load_phase(path_to_phase+str(i)+'.csv',lamb)
        psf_corono += apply_corono(phase_turb,psf_diff_limit,mask,step,pad)
        psf += compute_psf(phase_turb,mask,step,pad)
    psf_corono /= np.max(psf)



    return psf_corono

def plot_int_dd_comp(img_int, img_dd, D, n_pixels, lamb, step, pad):
    psf_size = int(n_pixels / step + 2 * pad)
    fov = n_pixels / step
    fs = fov
    res = fov/psf_size

    fov = n_pixels/step

    axis = np.arange(-fov / 2, fov / 2, res)
    labels = np.arange(-fov/2,fov/2+1, 10, dtype=np.int16)
    ticks = np.interp(labels,axis,np.arange(0,psf_size))

    vmin = 1e-6
    vmax = 3e-3

    plt.figure(figsize=(19, 7))
    plt.subplot(121)
    plt.imshow(img_int, norm=colors.LogNorm(vmin=vmin,vmax=vmax))
    plt.xticks(ticks, labels,fontsize=18)
    plt.yticks(ticks, labels,fontsize=18)
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label(label="log10 contrast", size=18)
    cbar.ax.tick_params(labelsize=18)
    plt.magma()
    plt.xlabel('Angular sepration [λ/D]', fontsize=18)
    plt.ylabel('Angular sepration [λ/D]', fontsize=18)
    plt.title('Integrator', fontsize=18)

    plt.subplot(122)
    plt.imshow(img_dd, norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    plt.xticks(ticks, labels,fontsize=18)
    plt.yticks(ticks, labels,fontsize=18)
    cbar = plt.colorbar(fraction=0.046, pad=0.04)
    cbar.set_label(label="log10 contrast", size=18)
    cbar.ax.tick_params(labelsize=18)
    plt.magma()
    plt.xlabel('Angular sepration [λ/D]', fontsize=18)
    plt.ylabel('Angular sepration [λ/D]', fontsize=18)
    plt.title('Data-driven', fontsize=18)

    plt.suptitle('Coronographic PSF comparison, wind = 34m/s', fontsize=18,fontweight="bold")
    plt.tight_layout(pad=0.0001)

    range = np.abs(axis)<30

    plt.figure(figsize=(19, 14))

    plt.subplot(211)
    cut_int = img_int[int(psf_size/2),:]
    cut_dd = img_dd[int(psf_size/2),:]

    plt.plot(axis[range],cut_int[range],axis[range],cut_dd[range])
    plt.yscale('log')
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize=18)
    plt.xlabel('Angular sepration [λ/D]', fontsize=18)
    plt.ylabel('log10 contrast', fontsize=18)
    plt.legend(['int','dd'],fontsize=18)
    plt.title('Cut parallel to wind axis', fontsize=18)
    plt.grid("True", which="both")

    plt.subplot(212)
    cut_int = img_int[:,int(psf_size/2)]
    cut_dd = img_dd[:,int(psf_size/2)]

    plt.plot(axis[range],cut_int[range],axis[range],cut_dd[range])
    plt.yscale('log')
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize=18)
    plt.xlabel('Angular sepration [λ/D]', fontsize=18)
    plt.ylabel('log10 contrast', fontsize=18)
    plt.legend(['int','dd'],fontsize=18)
    plt.title('Cut perpendicular to wind axis', fontsize=18)
    plt.grid("True",which="both")
    plt.tight_layout(pad = 6)
    plt.suptitle('PSF cut comparison, wind = 34m/s', fontsize=18, fontweight="bold")
    plt.show()
# def polt_radial(img_int):
#     dummy = hcipy.metrics.radial_profile(img_int,1)
#     return dummy

if __name__ == '__main__':
    lamb = 1.65
    n_frames = 996
    n_pixels = 400
    c_obs = 0.14
    step = 4
    pad = 50

    psf_corono_int = average_corono('phase_int_34/phase_tar_int_', 'phase_diff_limit/phase_diff_limit_tar.csv', lamb, n_frames, n_pixels, c_obs, step, pad)
    psf_corono_dd = average_corono('phase_dd_34/phase_tar_dd_', 'phase_diff_limit/phase_diff_limit_tar.csv', lamb, n_frames, n_pixels, c_obs, step, pad)

    plot_int_dd_comp(psf_corono_int,psf_corono_dd,8,n_pixels,lamb,step,pad)
    # polt_radial(psf_corono_int)
    print("fone")