import numpy as np
from matplotlib import pyplot as plt

class phase_plot:
    def __init__(self, title):
        self.fig, self.ax = plt.subplots(constrained_layout=True)
        self.fig.suptitle(title)

    def plot(self, phase, subtitle = ""):
        self.ax.imshow(phase)
        self.ax.set_title(subtitle)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()