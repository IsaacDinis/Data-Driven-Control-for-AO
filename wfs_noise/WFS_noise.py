import numpy as np
from matplotlib import pyplot as plt
from scipy.io import savemat, loadmat

def plot_res(t,res1,res2,res3):
    plt.figure(figsize=(7, 7))
    plt.plot(t,res2)
    plt.plot(t, res1)
    plt.plot(t,res3)
    plt.legend(['res1', 'res2','res3'])
    plt.show()


if __name__ == '__main__':
    fs = 1000
    wfs_mode = loadmat('data/single_mode_dist.mat')['data'].squeeze()
    wfs_mode_5 = loadmat('data/single_mode_dist_10.mat')['data'].squeeze()
    phase_mode = loadmat('data/single_mode_dist_phase.mat')['data'].squeeze()

    t = np.arange(0,wfs_mode.shape[0]/fs,1/fs)
    plot_res(t,wfs_mode,wfs_mode_5, phase_mode)