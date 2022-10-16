import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

def display_wf():
    mask = pupil_mask(404,0.14)
    inf_mat = np.genfromtxt('dummy.csv', delimiter=",")
    inf_mat = inf_mat*mask
    contrast = np.abs(np.fft.fft2(inf_mat))
    fig, ((a,b)) = plt.subplots(2,1)
    fig.suptitle("2 Frame Delay")
    im = a.imshow((inf_mat))#, norm = colors.LogNorm(vmin = 50, vmax=10000))
    a.title.set_text("int")
    im = b.imshow((contrast), norm = colors.LogNorm(vmin = 1, vmax=1000))
    b.title.set_text("dd")
    plt.show()

def pupil_mask(diam,c_obs):

    center = (int(diam/2), int(diam/2))

    r1 = center[0]
    r2 = diam*c_obs/2
    Y, X = np.ogrid[:diam, :diam]

    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = (dist_from_center <= r1) * (dist_from_center >= r2)
    return mask

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.

def contrast():
    mask = pupil_mask(404, 0.14)
    phase_dd_tot = np.zeros_like(np.abs(np.fft.fft2(np.genfromtxt('/home/isaac/shesha/residuals/phase_dd/phase_dd_0.csv', delimiter=","))))
    phase_int_tot = np.zeros_like(phase_dd_tot)

    for i in range(1,100):
        phase_dd = np.genfromtxt('/home/isaac/shesha/residuals/phase_dd/phase_dd_'+str(i)+'.csv', delimiter=",")
        phase_dd = phase_dd * mask
        phase_dd = np.abs(np.fft.fft2(phase_dd))
        phase_dd_tot += phase_dd

        phase_int = np.genfromtxt('/home/isaac/shesha/residuals/phase_int/phase_int_'+str(i)+'.csv', delimiter=",")
        phase_int = phase_int * mask
        phase_int = np.abs(np.fft.fft2(phase_int))
        phase_int_tot += phase_int

    fig, ((a, b)) = plt.subplots(2, 1)
    fig.suptitle("2 Frame Delay")
    im = a.imshow((phase_int_tot),norm=colors.LogNorm(vmin=1, vmax=10000))
    a.title.set_text("int")
    im = b.imshow((phase_dd_tot), norm=colors.LogNorm(vmin=1, vmax=10000))
    b.title.set_text("dd")
    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    contrast()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
