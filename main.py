import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

def display_wf():
    mask = pupil_mask(404,0.14)
    inf_mat = np.genfromtxt('phase_dd/phase_dd_99.csv', delimiter=",")
    inf_mat_masked = inf_mat*mask

    # contrast = np.abs(np.fft.fft2(inf_mat))
    fig, ((a,b)) = plt.subplots(2,1)
    fig.suptitle("Phase")
    im = a.imshow((inf_mat)) #, norm = colors.LogNorm(vmin = 1, vmax=1000))
    a.title.set_text("phase")
    im = b.imshow((inf_mat_masked))#, norm = colors.LogNorm(vmin = 1, vmax=1000))
    b.title.set_text("phase with mask")
    # plt.show()

def pupil_mask(diam,c_obs):

    center = (int(diam/2), int(diam/2))

    r1 = center[0]
    r2 = diam*c_obs/2
    Y, X = np.ogrid[:diam, :diam]

    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = 1.*((dist_from_center <= r1) * (dist_from_center >= r2))
    return mask

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.

def phase_psd(phase,mask,step,pad):
    phase -= np.mean(phase[mask > 0])
    phase = phase * mask
    phase = phase[::step, ::step]
    phase = np.pad(phase, pad, 'constant', constant_values=0)
    return np.fft.fftshift(np.abs(np.fft.fft2(phase)))

def intensity_coro(phase,strehl,mask,step,pad):
    phase /= 0.7/1.65
    phase -= np.mean(phase[mask > 0])

    e_field = np.exp(phase*1j)*mask

    e_field -= mask*np.sqrt(strehl)
    e_field = e_field[::step, ::step]
    e_field = np.pad(e_field, pad, 'constant', constant_values=0)

    return np.fft.fftshift(np.abs(np.fft.fft2(e_field))**2)

def psf(phase,mask,step,pad):
    phase /= 0.7/1.65
    phase -= np.mean(phase[mask > 0])

    e_field = np.exp(phase*1j)*mask

    e_field = e_field[::step, ::step]
    e_field = np.pad(e_field, pad, 'constant', constant_values=0)

    return np.fft.fftshift(np.abs(np.fft.fft2(e_field))**2)

def contrast():
    step = 8
    pad = 50
    mask = pupil_mask(404, 0.14)
    dummy = np.genfromtxt('phase_dd/phase_dd_3.csv', delimiter=",")
    dummy = np.pad(dummy, 400, 'constant', constant_values=0)
    phase_dd_tot = np.zeros_like(np.abs(np.fft.fft2(dummy[::step,::step])))
    phase_int_tot = np.zeros_like(phase_dd_tot)

    psf_int_tot = np.zeros_like(phase_dd_tot)
    psf_dd_tot = np.zeros_like(phase_dd_tot)

    strehl_int = np.genfromtxt('phase_int/strehl_int.csv', delimiter=",")
    strehl_dd = np.genfromtxt('phase_dd/strehl_dd.csv', delimiter=",")

    for i in range(3,100):
        phase_dd = np.genfromtxt('phase_dd/phase_dd_'+str(i)+'.csv', delimiter=",")
        # phase_dd = phase_psd(phase_dd, mask, step, 50)
        i_coro_dd = intensity_coro(phase_dd, strehl_dd[i], mask, step, pad)
        phase_dd_tot += i_coro_dd

        psf_dd = psf(phase_dd,mask,step,pad)
        psf_dd_tot += psf_dd

        phase_int = np.genfromtxt('phase_int/phase_int_'+str(i)+'.csv', delimiter=",")
        # phase_int = phase_psd(phase_int,mask,step,50)
        i_coro_int = intensity_coro(phase_int, strehl_int[i], mask, step, pad)
        phase_int_tot += i_coro_int

        psf_int = psf(phase_int,mask,step,pad)
        psf_int_tot += psf_int

    phase_dd_tot /= np.max(psf_dd_tot)
    phase_int_tot /= np.max(psf_int_tot)

    fig, ((a, b)) = plt.subplots(1, 2)
    fig.suptitle("phase PSD")
    im = a.imshow((phase_int_tot),norm=colors.LogNorm(vmin=2e-5, vmax=0.0028)) #250
    a.title.set_text("int")
    im = b.imshow((phase_dd_tot), norm=colors.LogNorm(vmin=2e-5, vmax=0.0028))
    b.title.set_text("dd")
    plt.show()

    # fig, ((a, b)) = plt.subplots(1, 2)
    # fig.suptitle("phase PSD")
    # im = a.imshow((phase_int_tot),norm=colors.LogNorm()) #(vmin=1, vmax=np.max(phase_int_tot))) #250
    # a.title.set_text("int")
    # im = b.imshow((phase_dd_tot), norm=colors.LogNorm())#min=1, vmax=np.max(phase_int_tot)))
    # b.title.set_text("dd")
    # plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    # display_wf()
    contrast()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
