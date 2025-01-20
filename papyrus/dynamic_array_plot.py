# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation


# class SingleArrayDynamicPlot:
#     def __init__(self, data_source, title="Dynamic Data Plot", x_label="Index", y_label="Value", interval=500,
#                  x_range=10):
#         """
#         Initialize the dynamic plot for a single array.

#         Parameters:
#         - data_source: object with a `get_data()` method that returns the array to be plotted
#         - title: title of the plot
#         - x_label: label for the x-axis
#         - y_label: label for the y-axis
#         - interval: update interval in milliseconds
#         - x_range: initial x-axis range
#         """
#         self.data_source = data_source
#         self.interval = interval
#         self.fig, self.ax = plt.subplots()
#         self.line, = self.ax.plot([], [], label=title)

#         self.ax.set_title(title)
#         self.ax.set_xlabel(x_label)
#         self.ax.set_ylabel(y_label)
#         self.ax.legend(loc='upper right')  # Set legend to North-East
#         self.ax.set_xlim(0, x_range)  # Set initial x-axis limits
#         self.init_ylim = (-1, 1)  # Default y-axis limits
#         self.ax.set_ylim(*self.init_ylim)

#     def init_plot(self):
#         """Initialize the plot."""
#         self.line.set_data([], [])
#         return self.line,

#     def update_plot(self, frame):
#         """Update the plot with new data."""
#         data_to_plot = self.data_source.get_data()

#         if len(data_to_plot) > 0:
#             # Update the line data
#             self.line.set_data(np.arange(len(data_to_plot)), data_to_plot)

#             # Dynamically update the y-axis limits
#             y_min, y_max = np.min(data_to_plot), np.max(data_to_plot)
#             self.ax.set_ylim(y_min - 0.1, y_max + 0.1)  # Add a small buffer

#             # Dynamically update the x-axis range if data grows
#             self.ax.set_xlim(0, len(data_to_plot))

#         return self.line,

#     def start(self):
#         """Start the animation."""
#         plt.ion()
#         self.ani = FuncAnimation(
#             self.fig,
#             self.update_plot,
#             init_func=self.init_plot,
#             blit=False,
#             interval=self.interval,
#             cache_frame_data=False
#         )
#         self.fig.show()


# class MultiArrayDynamicPlot:
#     def __init__(self, data_sources, labels, sampling_frequencies, x_range=10, title="Dynamic Data Plot",
#                  x_label="Time (s)", y_label="Value", interval=500):
#         """
#         Initialize the dynamic plot for multiple arrays with different sampling frequencies.

#         Parameters:
#         - data_sources: list of callables or objects, each with a `get_data()` method returning an array
#         - labels: list of strings, labels for each data series
#         - sampling_frequencies: list of floats, sampling frequencies (Hz) for each data source
#         - x_range: visible time range on the x-axis (in seconds)
#         - title: title of the plot
#         - x_label: label for the x-axis
#         - y_label: label for the y-axis
#         - interval: update interval in milliseconds
#         """
#         self.data_sources = data_sources
#         self.labels = labels
#         self.sampling_frequencies = sampling_frequencies
#         self.x_range = x_range
#         self.interval = interval
#         self.fig, self.ax = plt.subplots()
#         self.lines = [self.ax.plot([], [], label=label)[0] for label in labels]

#         self.ax.set_title(title)
#         self.ax.set_xlabel(x_label)
#         self.ax.set_ylabel(y_label)
#         self.ax.legend(loc='upper right')  # Set legend to North-East
#         self.ax.set_xlim(0, x_range)  # Enforce the x-axis range
#         self.init_ylim = (-1, 1)  # Default y-axis limits
#         self.ax.set_ylim(*self.init_ylim)

#     def resample_data(self, data, fs):
#         """
#         Resample the data to fit within the x_range.

#         Parameters:
#         - data: array-like, the data to be resampled
#         - fs: float, sampling frequency of the data

#         Returns:
#         - resampled_time: array, resampled time axis
#         - resampled_data: array, resampled data
#         """
#         num_points = int(self.x_range * fs)  # Number of points in x_range
#         if len(data) < num_points:
#             # Pad with zeros if data is shorter than required points
#             data = np.pad(data, (0, num_points - len(data)), 'constant')
#         else:
#             # Trim data to the last `num_points`
#             data = data[-num_points:]
#         time_axis = np.linspace(0, self.x_range, num_points)  # Uniformly distributed time axis
#         return time_axis, data

#     def init_plot(self):
#         """Initialize the plot."""
#         for line in self.lines:
#             line.set_data([], [])
#         return self.lines

#     def update_plot(self, frame):
#         """Update the plot with new data."""
#         all_y_data = []

#         for line, source, fs in zip(self.lines, self.data_sources, self.sampling_frequencies):
#             data_to_plot = source.get_data()
#             if len(data_to_plot) > 0:
#                 # Resample the data to fit within x_range
#                 time_axis, resampled_data = self.resample_data(data_to_plot, fs)

#                 # Update the line data
#                 line.set_data(time_axis, resampled_data)
#                 all_y_data.extend(resampled_data)

#         if all_y_data:
#             # Dynamically update the y-axis limits
#             y_min, y_max = np.min(all_y_data), np.max(all_y_data)
#             self.ax.set_ylim(y_min - 0.1, y_max + 0.1)  # Add a small buffer

#         return self.lines

#     def start(self):
#         """Start the animation."""
#         plt.ion()
#         self.ani = FuncAnimation(
#             self.fig,
#             self.update_plot,
#             init_func=self.init_plot,
#             blit=False,
#             interval=self.interval,
#             cache_frame_data=False
#         )
#         self.fig.show()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class SingleArrayDynamicPlot:
    def __init__(self, data_source, title="Dynamic Data Plot", x_label="Index", y_label="Value", interval=500,
                 x_range=10, mode_index=0,x_scale=None, y_scale=None):
        """
        Initialize the dynamic plot for a single array or a single mode from multi-dimensional data.

        Parameters:
        - data_source: object with a `get_data()` method that returns an array or a 2D array (x, n_modes)
        - title: title of the plot
        - x_label: label for the x-axis
        - y_label: label for the y-axis
        - interval: update interval in milliseconds
        - x_range: initial x-axis range
        - mode_index: integer, the mode to select if data is 2D (x, n_modes)
        """
        self.data_source = data_source
        self.interval = interval
        self.mode_index = mode_index
        self.fig, self.ax = plt.subplots()
        self.line, = self.ax.plot([], [], label=title)

        self.ax.set_title(title)
        self.ax.set_xlabel(x_label)
        self.ax.set_ylabel(y_label)
        self.ax.legend(loc='upper right')  # Set legend to North-East
        self.ax.set_xlim(0, x_range)  # Set initial x-axis limits
        self.init_ylim = (-1, 1)  # Default y-axis limits
        self.ax.set_ylim(*self.init_ylim)
        self.x_scale = x_scale
        self.y_scale = y_scale

    def init_plot(self):
        """Initialize the plot."""
        self.line.set_data([], [])
        return self.line,

    def update_plot(self, frame):
        """Update the plot with new data."""
        data_to_plot = self.data_source.get_data()

        if len(data_to_plot) > 0:
            # Handle multi-dimensional data by selecting the specified mode
            if data_to_plot.ndim == 2:
                data_to_plot = data_to_plot[:, self.mode_index]

            # Update the line data
            self.line.set_data(np.arange(len(data_to_plot)), data_to_plot)

            # Dynamically update the y-axis limits
            y_min, y_max = np.min(data_to_plot), np.max(data_to_plot)
            self.ax.set_ylim(y_min - 0.1, y_max + 0.1)  # Add a small buffer

            # Dynamically update the x-axis range if data grows
            self.ax.set_xlim(0, len(data_to_plot))
            
        if self.x_scale == "log":
            self.ax.set_xscale("log")
        if self.y_scale == "log":
            self.ax.set_yscale("log")

        return self.line,

    def start(self):
        """Start the animation."""
        plt.ion()
        self.ani = FuncAnimation(
            self.fig,
            self.update_plot,
            init_func=self.init_plot,
            blit=False,
            interval=self.interval,
            cache_frame_data=False
        )
        self.fig.show()


# class MultiArrayDynamicPlot:
#     def __init__(self, data_sources, labels, sampling_frequencies, x_range=10, title="Dynamic Data Plot",
#                  x_label="Time (s)", y_label="Value", interval=500, mode_indices=None,x_scale=None, y_scale=None):
#         """
#         Initialize the dynamic plot for multiple arrays with different sampling frequencies.

#         Parameters:
#         - data_sources: list of callables or objects, each with a `get_data()` method returning an array or 2D array
#         - labels: list of strings, labels for each data series
#         - sampling_frequencies: list of floats, sampling frequencies (Hz) for each data source
#         - x_range: visible time range on the x-axis (in seconds)
#         - title: title of the plot
#         - x_label: label for the x-axis
#         - y_label: label for the y-axis
#         - interval: update interval in milliseconds
#         - mode_indices: list of integers, modes to select for each data source if data is 2D
#         """
#         self.data_sources = data_sources
#         self.labels = labels
#         self.sampling_frequencies = sampling_frequencies
#         self.x_range = x_range
#         self.interval = interval
#         self.mode_indices = mode_indices or [0] * len(data_sources)  # Default to first mode for all sources
#         self.fig, self.ax = plt.subplots()
#         self.lines = [self.ax.plot([], [], label=label)[0] for label in labels]

#         self.ax.set_title(title)
#         self.ax.set_xlabel(x_label)
#         self.ax.set_ylabel(y_label)
#         self.ax.legend(loc='upper right')  # Set legend to North-East
#         self.ax.set_xlim(0, x_range)  # Enforce the x-axis range
#         self.init_ylim = (-1, 1)  # Default y-axis limits
#         self.ax.set_ylim(*self.init_ylim)
#         self.x_scale = x_scale
#         self.y_scale = y_scale

#     def resample_data(self, data, fs):
#         """
#         Resample the data to fit within the x_range.

#         Parameters:
#         - data: array-like, the data to be resampled
#         - fs: float, sampling frequency of the data

#         Returns:
#         - resampled_time: array, resampled time axis
#         - resampled_data: array, resampled data
#         """
#         num_points = int(self.x_range * fs)  # Number of points in x_range
#         if len(data) < num_points:
#             # Pad with zeros if data is shorter than required points
#             data = np.pad(data, (0, num_points - len(data)), 'constant')
#         else:
#             # Trim data to the last `num_points`
#             data = data[-num_points:]
#         time_axis = np.linspace(0, self.x_range, num_points)  # Uniformly distributed time axis
#         return time_axis, data

#     def init_plot(self):
#         """Initialize the plot."""
#         for line in self.lines:
#             line.set_data([], [])
#         return self.lines

#     def update_plot(self, frame):
#         """Update the plot with new data."""
#         all_y_data = []

#         for line, source, fs, mode_index in zip(self.lines, self.data_sources, self.sampling_frequencies, self.mode_indices):
#             data_to_plot = source.get_data()
#             if len(data_to_plot) > 0:
#                 # Handle multi-dimensional data by selecting the specified mode
#                 if data_to_plot.ndim == 2:
#                     data_to_plot = data_to_plot[:, mode_index]

#                 # Resample the data to fit within x_range
#                 time_axis, resampled_data = self.resample_data(data_to_plot, fs)

#                 # Update the line data
#                 line.set_data(time_axis, resampled_data)
#                 all_y_data.extend(resampled_data)

#         if all_y_data:
#             # Dynamically update the y-axis limits
#             y_min, y_max = np.min(all_y_data), np.max(all_y_data)
#             self.ax.set_ylim(y_min - 0.1, y_max + 0.1)  # Add a small buffer
#         if self.y_scale == "log":
#             if np.any(np.array(all_y_data) <= 0):
#                 raise ValueError("Log scale requires all y-values to be strictly positive.")
#             self.ax.set_yscale("log", nonpositive='clip')

#         if self.x_scale == "log":
#             self.ax.set_xscale("log", nonpositive='clip')

#         return self.lines

#     def start(self):
#         """Start the animation."""
#         plt.ion()
#         self.ani = FuncAnimation(
#             self.fig,
#             self.update_plot,
#             init_func=self.init_plot,
#             blit=False,
#             interval=self.interval,
#             cache_frame_data=False
#         )
#         self.fig.show()

class MultiArrayDynamicPlot:
    def __init__(self, data_sources, labels, x_sources, title="Dynamic Data Plot",
                 x_label="Time (s)", y_label="Value", interval=500, mode_indices=None, x_scale=None, y_scale=None):
        """
        Initialize the dynamic plot for multiple arrays.

        Parameters:
        - data_sources: list of callables or objects, each with a `get_data()` method returning an array or 2D array
        - labels: list of strings, labels for each data series
        - x_sources: list of callables or objects, each with a `get_x_data()` method returning the x-axis values
        - title: title of the plot
        - x_label: label for the x-axis
        - y_label: label for the y-axis
        - interval: update interval in milliseconds
        - mode_indices: list of integers, modes to select for each data source if data is 2D
        - x_scale: string, scale of the x-axis ("linear" or "log")
        - y_scale: string, scale of the y-axis ("linear" or "log")
        """
        self.data_sources = data_sources
        self.labels = labels
        self.x_sources = x_sources
        self.interval = interval
        self.mode_indices = mode_indices or [0] * len(data_sources)  # Default to first mode for all sources
        self.fig, self.ax = plt.subplots()
        self.lines = [self.ax.plot([], [], label=label)[0] for label in labels]

        self.ax.set_title(title)
        self.ax.set_xlabel(x_label)
        self.ax.set_ylabel(y_label)
        self.ax.legend(loc='upper right')
        self.init_ylim = (-1, 1)  # Default y-axis limits
        self.ax.set_ylim(*self.init_ylim)
        self.x_scale = x_scale
        self.y_scale = y_scale

    def init_plot(self):
        """Initialize the plot."""
        for line in self.lines:
            line.set_data([], [])
        return self.lines

    def update_plot(self, frame):
        """Update the plot with new data."""
        all_y_data = []
        all_x_data = []

        for line, source, x_source, mode_index in zip(self.lines, self.data_sources, self.x_sources, self.mode_indices):
            data_to_plot = source.get_data()
            x_axis = x_source.get_data()
            # print(len(data_to_plot))
            if len(data_to_plot) > 0 and len(x_axis) > 0:
                # Handle multi-dimensional data by selecting the specified mode
                if data_to_plot.ndim == 2:
                    data_to_plot = data_to_plot[:, mode_index]

                # Ensure data for log-scale compliance
                # if self.y_scale == "log":
                #     data_to_plot = np.clip(data_to_plot, 1e-10, None)
                # if self.x_scale == "log":
                #     x_axis = np.clip(x_axis, 1e-10, None)

                # Update the line data
                line.set_data(x_axis, data_to_plot)
                all_y_data.extend(data_to_plot)
                all_x_data.extend(x_axis)

        if all_y_data and all_x_data:
            # Dynamically adjust y-axis limits
            if self.y_scale == "log":
                all_y_data = np.clip(all_y_data,1e-5,np.inf)
            y_min, y_max = np.min(all_y_data), np.max(all_y_data)
            # if y_min <= 0 and self.y_scale == "log":
                # y_min = max(y_min, 1e-3)  # Ensure positive y_min for log scale

            # else:
            #     self.ax.set_ylim(y_min - 0.1, y_max + 0.1)  # Add a small buffer

            # Dynamically adjust x-axis limits
            x_min, x_max = np.min(all_x_data), np.max(all_x_data)
            self.ax.set_xlim(x_min, x_max)  # Add a small buffer
            self.ax.set_ylim(y_min, y_max)

        if self.x_scale == "log":
            self.ax.set_xscale("log")

        if self.y_scale == "log":
            
            self.ax.set_yscale("log")
            print(np.min(all_y_data))
            # print(len(all_y_data))
            # print(len(all_x_data))

        return self.lines

    def start(self):
        """Start the animation."""
        plt.ion()
        self.ani = FuncAnimation(
            self.fig,
            self.update_plot,
            init_func=self.init_plot,
            blit=False,
            interval=self.interval,
            cache_frame_data=False
        )
        self.fig.show()