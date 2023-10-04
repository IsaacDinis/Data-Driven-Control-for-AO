import numpy as np
from matplotlib import pyplot as plt
from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 


class zernike_plot:
    def __init__(self, title, refresh_rate, order, pup_diam):
        pupil_grid = make_pupil_grid(pup_diam)
        zernike_basis = make_zernike_basis(order, 1, pupil_grid)
        transform_matrix = zernike_basis.transformation_matrix
        transform_matrix = np.linalg.pinv(transform_matrix)
        # self.order = order
        self.transform_matrix = transform_matrix[1:,:]
        self.fig, self.ax = plt.subplots(constrained_layout=True)
        self.fig.suptitle(title)
        self.refresh_rate = refresh_rate
        self.count = 0
        self.zernike_res = np.zeros(order-1)
        # self.line, = self.ax.plot(self.zernike_res)
        self.ax.stem(self.zernike_res)
        self.ax.set_ylabel("[nm]")
        self.ax.set_xlabel("order")
        self.dummy = 0

    def plot(self, phase,iter_n):
        iter_n += 1
        if np.sum([phase!=0]) !=0:
            self.zernike_res += np.abs(self.transform_matrix@phase.flatten())
            self.dummy += np.abs(phase.flatten()/np.sum([phase!=0]))
        # self.zernike_res += np.square(phase.flatten()/np.sum([phase!=0]))
        zernike_rms = self.zernike_res/iter_n*1e3
        subtitle = 'tot = {:.3f}'.format(np.sum(zernike_rms))
        # subtitle = 'tot = {:.3f}'.format(np.sum(self.dummy)/iter_n*1e3)
        if self.count == 0:
            # self.line.set_ydata(zernike_rms)
            self.ax.cla()
            self.ax.stem(zernike_rms)
            self.ax.set_ylim(0,np.max(zernike_rms))
            # self.ax.autoscale()
            self.ax.set_title(subtitle)
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        self.count += 1
        if self.count == self.refresh_rate:
            self.count = 0

    def reset(self):
        self.zernike_res = 0
        self.dummy = 0