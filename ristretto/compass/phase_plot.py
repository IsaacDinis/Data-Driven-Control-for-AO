import numpy as np
from matplotlib import pyplot as plt

class phase_plot:
    def __init__(self):
        self.fig, self.ax = plt.subplots()
    def plot(self, phase):
        self.ax.imshow(phase)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()