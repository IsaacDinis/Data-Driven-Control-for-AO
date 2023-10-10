import numpy as np
from matplotlib import pyplot as plt
import compute_psd
import astropy.io.fits as pfits

class modal_plot:
    def __init__(self, title, refresh_rate, order, n_iter):

        self.fig, self.ax = plt.subplots(constrained_layout=True)
        self.fig.suptitle(title)
        self.refresh_rate = refresh_rate
        self.count = 0
        # self.line, = self.ax.plot(self.modal_res)
        self.ax.stem(np.zeros(order))
        self.ax.set_ylabel("[au]")
        self.ax.set_xlabel("order")
        self.modal_res = np.zeros((n_iter,order))
        self.n_iter = n_iter

    def plot(self, modal_res,iter_n, bool_update = True):
        if bool_update:
            self.modal_res[iter_n,:] = modal_res

        subtitle = 'tot = {:.3f}'.format(np.sqrt(np.sum(np.square(modal_res))))

        if self.count == 0:
            # self.line.set_ydata(modal_rms)
            self.ax.cla()
            self.ax.stem(np.abs(modal_res),'-')
            self.ax.set_ylim(0,np.max(np.abs(modal_res)))
            self.ax.set_title(subtitle)
            self.ax.set_ylabel("[au]")
            self.ax.set_xlabel("order")
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        self.count += 1
        if self.count == self.refresh_rate:
            self.count = 0

    def reset(self):
        self.modal_res *= 0

    def plot_psd(self,mode,fs,title,path_name = ''):
        n_average = 10
        window_size = int(np.floor(self.n_iter/n_average))
        psd, freq = compute_psd.compute_psd_fft(self.modal_res, n_average, window_size,fs)
        psd /= np.max(psd)
        psd = np.log10(psd)
        fig, ax = plt.subplots(constrained_layout=True)
        fig.suptitle(title)
        ax.set_title('mode {:d}'.format(mode))
        ax.set_xscale('log')
        ax.set_ylabel("mag. [dB]")
        ax.set_xlabel("freq. [Hz]")
        ax.plot(freq,psd[:,mode])
        if (path_name != ''):
            fig.savefig(path_name) 

    def save(self, path_name):
        pfits.writeto(path_name, self.modal_res, overwrite = True)

    def load(self,path_name):
        self.modal_res = pfits.getdata(path_name)

    def save_std_plot(self, path_name):
        modal_res_rms = np.std(self.modal_res,axis = 0)
        self.plot(modal_res_rms,0,False)
        self.fig.savefig(path_name) 