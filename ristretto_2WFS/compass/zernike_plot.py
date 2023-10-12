import numpy as np
from matplotlib import pyplot as plt
from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 
import compute_psd
import astropy.io.fits as pfits

class zernike_plot:
    def __init__(self, title, refresh_rate, order, pup_diam, pupil, n_iter):
        pupil_grid = make_pupil_grid(pup_diam)
        zernike_basis = make_zernike_basis(order, 1, pupil_grid)
        transform_matrix = zernike_basis.transformation_matrix
        transform_matrix = np.linalg.pinv(transform_matrix)
        transform_matrix /= np.std(transform_matrix[1,:],where = transform_matrix[1,:]!=0)
        transform_matrix /= np.sum(pupil)
        # self.order = order
        self.transform_matrix = transform_matrix[1:,:]
        self.fig, self.ax = plt.subplots(constrained_layout=True)
        self.fig.suptitle(title)
        self.refresh_rate = refresh_rate
        self.count = 0
        
        # self.line, = self.ax.plot(self.zernike_res)
        self.ax.stem(np.zeros(order-1))
        self.ax.set_ylabel("[nm]")
        self.ax.set_xlabel("order")
        self.zernike_res = np.zeros((n_iter,order-1))
        self.n_iter = n_iter


    def plot(self, phase,iter_n):
        
        if np.sum([phase!=0]) !=0:
            zernike_res = self.transform_matrix@phase.flatten()*1e3
        else:
            zernike_res = np.zeros_like(self.zernike_res[0,:])
        self.zernike_res[iter_n,:] = zernike_res

        subtitle = 'tot = {:.3f}'.format((np.sqrt(np.sum(np.square(zernike_res)))))

        if self.count == 0:
            # self.line.set_ydata(zernike_rms)
            self.ax.cla()
            self.ax.stem(np.abs(zernike_res))
            self.ax.set_ylim(0,np.max(np.abs(zernike_res)))
            self.ax.set_title(subtitle)
            self.ax.set_ylabel("[nm]")
            self.ax.set_xlabel("order")
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        self.count += 1
        if self.count == self.refresh_rate:
            self.count = 0

    def reset(self):
        self.zernike_res *= 0


    def plot_psd(self,mode,fs,title,path_name = ''):
        n_average = 10
        window_size = int(np.floor(self.n_iter/n_average))
        psd, freq = compute_psd.compute_psd_fft(self.zernike_res, n_average, window_size,fs)
        fig, ax = plt.subplots(constrained_layout=True)
        # psd /= np.max(psd)
        psd = np.log10(psd)
        fig.suptitle(title)
        ax.set_title('mode {:d}'.format(mode))
        ax.set_xscale('log')
        ax.set_ylabel("mag. [dB]")
        ax.set_xlabel("freq. [Hz]")
        ax.plot(freq,psd[:,mode])
        if (path_name != ''):
            fig.savefig(path_name) 

    def save(self, path_name):
        pfits.writeto(path_name, self.zernike_res, overwrite = True)

    def load(self,path_name):
        self.zernike_res = pfits.getdata(path_name)
        
    def save_std_plot(self, path_name):
        zernike_res_rms = np.std(self.zernike_res,axis = 0)
        subtitle = 'tot = {:.3f}'.format((np.sqrt(np.sum(np.square(zernike_res_rms)))))
        self.ax.cla()
        self.ax.stem(zernike_res_rms)
        self.ax.set_ylim(0,np.max(zernike_res_rms))
        self.ax.set_title(subtitle)
        self.ax.set_ylabel("[nm]")
        self.ax.set_xlabel("order")
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        self.fig.savefig(path_name) 
