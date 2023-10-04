import numpy as np
from matplotlib import pyplot as plt

class phase_plot:
    def __init__(self, title, refresh_rate):
        self.fig, self.ax = plt.subplots(constrained_layout=True)
        self.im = self.ax.imshow(np.zeros((8,8)))
        self.cbar = self.fig.colorbar(self.im, ax=self.ax)
        self.cbar.set_label(label="[um]", size=12)
        self.fig.suptitle(title)
        self.refresh_rate = refresh_rate
        self.count = 0
        self.opd = 0
        self.opd_rms = 0

    def plot(self, phase, subtitle = "", iter_n = 0):
        # self.opd += np.sum(np.square(phase)/np.sum([phase!=0]))
        if np.sum([phase!=0]) !=0:
            self.opd += np.sum(np.abs(phase)/np.sum([phase!=0]))
        iter_n += 1
        # self.opd_rms = np.sqrt(self.opd/iter_n)*1e3
        self.opd_rms = (self.opd/iter_n)*1e3

        if self.count == 0:
            # im = self.ax.imshow(phase)
            self.im.set_data(phase)
            self.im.autoscale()

            self.ax.set_title(subtitle)
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        self.count += 1
        if self.count == self.refresh_rate:
            self.count = 0

    def reset(self):
        self.opd = 0
        self.opd_rms = 0
