import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

def display_wf():
    mask = pupil_mask(404,0.14)
    inf_mat = np.genfromtxt('phase_int/phase_int_3.csv', delimiter=",")
    inf_mat_masked = inf_mat*mask

    # contrast = np.abs(np.fft.fft2(inf_mat))
    fig, ((a,b)) = plt.subplots(2,1)
    fig.suptitle("Phase")
    im = a.imshow((inf_mat)) #, norm = colors.LogNorm(vmin = 1, vmax=1000))
    a.title.set_text("phase")
    # plt.colorbar(im, cax=a, orientation='vertical')
    # im = b.imshow((inf_mat_masked))#, norm = colors.LogNorm(vmin = 1, vmax=1000))
    im = b.imshow((mask>0))  #
    b.title.set_text("phase with mask")
    plt.show()

def pupil_mask(pupil_diam, c_obs):

    center = (int(pupil_diam / 2), int(pupil_diam / 2))

    r1 = center[0]
    r2 = pupil_diam * c_obs / 2
    Y, X = np.ogrid[:pupil_diam, :pupil_diam]

    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = 1.*((dist_from_center <= r1) * (dist_from_center >= r2))
    return mask

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.
def compute_strehl(phase_turb,phase_diff_limit,c_obs):
    mask = pupil_mask(phase_turb.shape[0],c_obs)
    psf_turb = compute_psf(phase_turb,mask)
    psf_diff_limit = compute_psf(phase_diff_limit,mask)
    strehl = np.max(psf_turb)/np.max(psf_diff_limit)
    return strehl

def phase_psd(phase,mask,step,pad):
    phase -= np.mean(phase[mask > 0])
    phase = phase * mask
    phase = phase[::step, ::step]
    phase = np.pad(phase, pad, 'constant', constant_values=0)
    return np.fft.fftshift(np.abs(np.fft.fft2(phase)))

def intensity_coro(phase,strehl,mask,step=1,pad=0,phase_factor=1,strehl_factor=1):
    phase *= phase_factor
    phase -= np.mean(phase[mask > 0])

    e_field = np.exp(phase*1j)*mask
    phase_diff_limit = np.genfromtxt('phase_diff_limit/phase_diff_limit_tar.csv', delimiter=",")
    strehl = compute_strehl(phase,phase_diff_limit,0.14)
    e_field -= mask*np.sqrt(strehl_factor*strehl)
    e_field = e_field[::step, ::step]
    e_field = np.pad(e_field, pad, 'constant', constant_values=0)

    return np.fft.fftshift(np.abs(np.fft.fft2(e_field))**2)

def compute_psf(phase, mask, step=1, pad=0):
    # phase *= phase_factor
    phase -= np.mean(phase[mask > 0])
    print(np.exp(-np.std(phase[mask>0])**2))
    e_field = np.exp(phase*1j)*mask


    e_field = e_field[::step, ::step]
    e_field = np.pad(e_field, pad, 'constant', constant_values=0)

    plt.imshow(np.angle(e_field))
    plt.colorbar()
    plt.show()


    psf = np.fft.fftshift(np.abs(np.fft.fft2(e_field))**2)
    return psf


def plot_int_dd_comp(img_int,img_dd,D,pupil_diam,lamb,step,pad):
    fs = pupil_diam/(step*D)
    res = fs/((pupil_diam/step)+2*pad)
    res /= lamb
    # N =
    # ticks =
    fig, ((a, b)) = plt.subplots(1, 2)
    fig.suptitle("phase PSD")
    a.imshow((img_int),norm=colors.LogNorm())#vmin=2e-5, vmax=0.0028)) #250
    a.title.set_text("int")
    ticks = np.arange(-fs/2,fs/2,res)/lamb
    a.set_xticks([1.5])
    a.set_xticklabels([7])
    # a.set_xticklabels(ticks)


    b.imshow((img_dd), norm=colors.LogNorm())#vmin=2e-5, vmax=0.0028))
    b.title.set_text("dd")
    # b.set_xticklabels(ticks)

def contrast():
    step = 8
    pad = 50
    mask_wfs = pupil_mask(404, 0.14)
    mask_tar = pupil_mask(400, 0.14)
    dummy = np.genfromtxt('phase_int/phase_int_67.csv', delimiter=",")
    dummy2 = np.genfromtxt('phase_int/tar_int_67.csv', delimiter=",")
    # plt.imshow(dummy*mask_wfs)
    # plt.colorbar()
    # plt.show()
    dummy3 = np.divide(dummy[2:-2,2:-2],dummy2)
    np.exp(-np.std(dummy[mask_wfs>0])**2)
    dummy = np.pad(dummy, 400, 'constant', constant_values=0)
    dummy2 = np.pad(dummy2, 400, 'constant', constant_values=0)

    phase_dd_tot = np.zeros_like(np.abs(np.fft.fft2(dummy[::step,::step])))
    phase_int_tot = np.zeros_like(phase_dd_tot)

    phase_dd_tar_tot = np.zeros_like(np.abs(np.fft.fft2(dummy2[::step, ::step])))
    phase_int_tar_tot = np.zeros_like(phase_dd_tar_tot)

    psf_int_tot = np.zeros_like(phase_dd_tot)
    psf_dd_tot = np.zeros_like(phase_dd_tot)

    strehl_int = np.genfromtxt('phase_int/strehl_int.csv', delimiter=",")
    strehl_dd = np.genfromtxt('phase_dd/strehl_dd.csv', delimiter=",")

    for i in range(3,100):
        # phase_dd = np.genfromtxt('phase_dd/phase_dd_'+str(i)+'.csv', delimiter=",")
        # # phase_dd = phase_psd(phase_dd, mask_wfs, step, 50)
        # i_coro_dd = intensity_coro(phase_dd, strehl_dd[i], mask_wfs, step, pad,1)
        # phase_dd_tot += i_coro_dd
        #
        # psf_dd = compute_psf(phase_dd, mask_wfs, step, pad)
        # psf_dd_tot += psf_dd
        #
        # phase_int = np.genfromtxt('phase_int/phase_int_'+str(i)+'.csv', delimiter=",")
        # # phase_int = phase_psd(phase_int,mask_wfs,step,50)
        # i_coro_int = intensity_coro(phase_int, strehl_int[i], mask_wfs, step, pad,1)
        # phase_int_tot += i_coro_int
        #
        # psf_int = compute_psf(phase_int, mask_wfs, step, pad)
        # psf_int_tot += psf_int


        phase_dd_tar = np.genfromtxt('phase_dd/tar_dd_'+str(i)+'.csv', delimiter=",")
        # phase_dd = phase_psd(phase_dd, mask_wfs, step, 50)
        i_coro_dd_tar = intensity_coro(phase_dd_tar, strehl_dd[i], mask_tar, step, pad,1,1)
        phase_dd_tar_tot += i_coro_dd_tar

        phase_int_tar = np.genfromtxt('phase_int/tar_int_'+str(i)+'.csv', delimiter=",")
        # phase_int = phase_psd(phase_int,mask_wfs,step,50)
        i_coro_int_tar = intensity_coro(phase_int_tar, strehl_int[i], mask_tar, step, pad,1)
        phase_int_tar_tot += i_coro_int_tar


    # phase_dd_tot /= np.max(psf_dd_tot)
    # phase_int_tot /= np.max(psf_int_tot)

    # plot_int_dd_comp(phase_int_tot, phase_dd_tot, 8, 404, 1, step, pad)
    plot_int_dd_comp(phase_int_tar_tot, phase_dd_tar_tot, 8, 400, 1, step, pad)
    # dummy3 = np.divide(phase_int_tot[1:,1:],phase_int_tar_tot)
    # fig, ((a, b)) = plt.subplots(1, 2)
    # fig.suptitle("phase PSD")
    # im = a.imshow((phase_int_tot),norm=colors.LogNorm()) #(vmin=1, vmax=np.max(phase_int_tot))) #250
    # a.title.set_text("int")
    # im = b.imshow((phase_dd_tot), norm=colors.LogNorm())#min=1, vmax=np.max(phase_int_tot)))
    # b.title.set_text("dd")
    plt.show()

def average_var():
    mask = pupil_mask(400,0.14)
    var = 0
    count = 0
    for i in range(3,100):
        count += 1
        phase = np.genfromtxt('phase_turb/phase_turb_tar_'+str(i)+'.csv', delimiter=",")
        phase *= 2*np.pi/(0.7)
        phase -= np.mean(phase[mask > 0])
        var += np.std(phase[mask > 0]) ** 2
        # phase_dd = phase_psd(phase_dd, mask_wfs, step, 50)
    var/= count
    print(var)
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    # display_wf()
    # average_var()

    mask = pupil_mask(400,0.14)

    # phase_turb = np.genfromtxt('phase_diff_limit/phase_turb_tar.csv', delimiter=",")
    phase_turb = np.genfromtxt('phase_diff_limit/phase_turb_tar.csv', delimiter=",")*2*np.pi/(0.7)
    phase_diff_limit = np.genfromtxt('phase_diff_limit/phase_diff_limit_tar.csv', delimiter=",")*2*np.pi/(0.7)
    psf_turb = compute_psf(phase_turb, mask)
    psf_diff_limit = compute_psf(phase_diff_limit, mask)
    print(np.max(psf_turb)/np.max(psf_diff_limit))
    # plt.imshow(psf_turb)#,norm=colors.LogNorm())
    # plt.colorbar()
    plt.plot(psf_turb[:,int(psf_turb.shape[0]/2)])
    # plt.yscale("log")
    plt.show()

    # contrast()
    # phase_diff_limit = np.genfromtxt('phase_diff_limit/phase_diff_limit_tar.csv', delimiter=",")
    # phase_turb = np.genfromtxt('phase_int/tar_int_67.csv', delimiter=",")
    # strehl = compute_strehl(phase_turb,phase_diff_limit,0.14)
    # print(np.std(phase_turb))
    # # print(strehl)
    #
    # phase_diff_limit = np.genfromtxt('phase_diff_limit/phase_diff_limit_wfs.csv', delimiter=",")
    # phase_turb = np.genfromtxt('phase_int/phase_int_67.csv', delimiter=",")
    # print(np.std(phase_turb))
    # strehl = compute_strehl(phase_turb,phase_diff_limit,0.14)

    # print(strehl)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
