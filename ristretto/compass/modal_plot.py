import numpy as np
from matplotlib import pyplot as plt



class modal_plot:
    def __init__(self, title, refresh_rate, order):

        self.fig, self.ax = plt.subplots(constrained_layout=True)
        self.fig.suptitle(title)
        self.refresh_rate = refresh_rate
        self.count = 0
        self.modal_res = np.zeros(order)
        # self.line, = self.ax.plot(self.modal_res)
        self.ax.stem(self.modal_res)
        self.ax.set_ylabel("[au]")
        self.ax.set_xlabel("order")


    def plot(self, modal_res,iter_n):
        iter_n += 1
        self.modal_res += np.square(modal_res)
        
        modal_rms = np.sqrt(self.modal_res/iter_n)
        subtitle = 'tot = {:.3f}'.format(np.sum(modal_rms))

        if self.count == 0:
            # self.line.set_ydata(modal_rms)
            self.ax.cla()
            self.ax.stem(modal_rms,'-')
            self.ax.set_ylim(0,np.max(modal_rms))
            self.ax.set_title(subtitle)
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        self.count += 1
        if self.count == self.refresh_rate:
            self.count = 0

    def reset(self):
        self.modal_res = 0
