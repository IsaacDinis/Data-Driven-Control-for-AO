from dynamic_array_plot import SingleArrayDynamicPlot, MultiArrayDynamicPlot
import dao

class SharedMemorySource:
    """A wrapper class for shared memory objects."""
    def __init__(self, shm_path):
        self.shm = dao.shm(shm_path)

    def get_data(self):
        """Fetch data from the shared memory."""
        return self.shm.get_data()[:, 0]

turb_source = SharedMemorySource('/tmp/turb_buf.shm')
pol_source = SharedMemorySource('/tmp/pol_buf.shm')
res_source = SharedMemorySource('/tmp/res_buf.shm')
command_source = SharedMemorySource('/tmp/command_buf.shm')

turb_fft_source = SharedMemorySource('/tmp/turb_buf_fft.shm')
pol_fft_source = SharedMemorySource('/tmp/pol_buf_fft.shm')
res_fft_source = SharedMemorySource('/tmp/res_buf_fft.shm')
command_fft_source = SharedMemorySource('/tmp/command_buf_fft.shm')
t_source = SharedMemorySource('/tmp/t.shm')
f_source = SharedMemorySource('/tmp/f_buff.shm')

# pol_plot = MultiArrayDynamicPlot(
#     data_sources=[turb_source, pol_source,],
#     labels=["turb", "pol"],
#     sampling_frequencies= [100,100],
#     title="pol",
#     x_label="time",
#     y_label="Value",
#     interval=1000,
#     x_range=10.24
# )

# ao_plot = MultiArrayDynamicPlot(
#     data_sources=[turb_source, res_source, command_source],
#     labels=["turb", "res", "command"],
#     sampling_frequencies= [100,100,100],
#     title="ao plot",
#     x_label="time",
#     y_label="Value",
#     interval=1000,
#     x_range=10.24
# )

# fft_plot = MultiArrayDynamicPlot(
#     data_sources=[turb_fft_source, res_fft_source, res_fft_source],
#     labels=["turb", "res", "command"],
#     sampling_frequencies= [1,1,1],
#     title="fft plot",
#     x_label="time",
#     y_label="Value",
#     interval=1000,
#     x_range=500,
#     y_scale = 'log')

# pol_plot.start()
# ao_plot.start()
# fft_plot.start()

    # x_sources = [DummyXSource()]
    # sources = [DummySource()]
    # labels = ["Sine Wave"]

    # plot = MultiArrayDynamicPlot(sources, labels, x_sources, x_scale="linear", y_scale="linear")
    # plot.start()

pol_plot = MultiArrayDynamicPlot(
    data_sources=[turb_source, pol_source,],
    labels=["turb", "pol"],
    x_sources = [t_source,t_source],
    title="pol",
    x_label="time",
    y_label="Value",
    interval=1000,
)

ao_plot = MultiArrayDynamicPlot(
    data_sources=[turb_source, res_source, command_source],
    labels=["turb", "res", "command"],
    x_sources = [t_source,t_source,t_source],
    title="ao plot",
    x_label="time",
    y_label="Value",
    interval=1000,
)

# fft_plot = MultiArrayDynamicPlot(
#     data_sources=[turb_fft_source, res_fft_source, res_fft_source],
#     labels=["turb", "res", "command"],
#     x_sources = [f_source,f_source,f_source],
#     title="fft plot",
#     x_label="time",
#     y_label="Value",
#     interval=1000,
#     y_scale = 'log')

fft_plot = MultiArrayDynamicPlot(
    data_sources=[pol_fft_source,res_fft_source,command_fft_source],
    labels=["pol",'res','command'],
    x_sources = [f_source,f_source,f_source],
    title="fft plot",
    x_label="time",
    y_label="Value",
    interval=1000,
    y_scale = 'log',
    x_scale = 'log')

pol_fft_plot = MultiArrayDynamicPlot(
    data_sources=[turb_fft_source,pol_fft_source],
    labels=["turb",'pol'],
    x_sources = [f_source,f_source],
    title="pol_fft plot",
    x_label="time",
    y_label="Value",
    interval=1000,
    y_scale = 'log',
    x_scale = 'log')


pol_plot.start()
fft_plot.start()
pol_fft_plot.start()
ao_plot.start()