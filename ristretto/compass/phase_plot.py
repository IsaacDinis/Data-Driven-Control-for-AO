import numpy as np
from matplotlib import pyplot as plt

class phase_plot:
    def __init__(self, title, refresh_rate):
        self.fig, self.ax = plt.subplots(constrained_layout=True)
        self.fig.suptitle(title)
        self.refresh_rate = refresh_rate
        self.count = 0

    def plot(self, phase, subtitle = ""):
        if self.count == 0:
            self.ax.imshow(phase)
            self.ax.set_title(subtitle)
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        self.count += 1
        if self.count == self.refresh_rate:
            self.count = 0