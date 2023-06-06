"""
script to control one mode (tilt)

Usage:
  closed_loop_tilt.py <parameters_filename> [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
  -i, --interactive  keep the script interactive
  -d, --devices devices      Specify the devices
  -n, --niter niter       Number of iterations
"""

from shesha.config import ParamConfig
from docopt import docopt
import numpy as np
from scipy.io import savemat, loadmat
import astropy.io.fits as pfits
import matplotlib.pyplot as plt
#ipython -i shesha/widgets/widget_ao.py ~/Data-Driven-Control-for-AO/2DM_study/compass/compass_param.py
#V2V = np.load('../../saxo-plus/Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/V_DM0_2_V_DM1.npy')

if __name__ == "__main__":
    arguments = docopt(__doc__)

    param_file = arguments["<parameters_filename>"]

    config = ParamConfig(param_file)

    from shesha.supervisor.compassSupervisor import CompassSupervisor as Supervisor

    if arguments["--devices"]:
        config.p_loop.set_devices([
                int(device) for device in arguments["--devices"].split(",")
        ])

    if arguments["--niter"]:
        n_iter = (int(arguments["--niter"]))
    else:
        n_iter = 4000


    supervisor = Supervisor(config)
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(False) 

    n_modes_DM0 = 88
    n_modes_DM1 = 800

    a = np.array([1.,-1]) 
    b = np.array([0.5,0])


    # Load command and influence matrix
    S2M_DM0 = np.load('calib_mat/S2M_DM0.npy')
    S2M_DM1 = np.load('calib_mat/S2M_DM1.npy')

    M2V_DM0 = np.load('calib_mat/M2V_DM0.npy')
    M2V_DM1 = np.load('calib_mat/M2V_DM1.npy')

    V_DM0_2_V_DM1 = np.load('calib_mat/V_DM0_2_V_DM1.npy')
    M_DM0_2_M_DM1 = np.load('calib_mat/M_DM0_2_M_DM1.npy')

    n_points = 1000
    res_DM0_0 = np.zeros(n_points)
    res_DM1_0 = np.zeros(n_points)

    res_DM0_1 = np.zeros(n_points)
    res_DM1_1 = np.zeros(n_points)

    #------------------------------------
    # control tilt mode
    #------------------------------------


    amp = np.linspace(-10,10,n_points)

    for i in range(n_points):

        voltage_DM0 = M2V_DM0[:,0]*amp[i]
        voltage = np.concatenate((voltage_DM0, np.zeros(M2V_DM1.shape[0])), axis=0)
        supervisor.rtc.set_command(0, voltage)

        supervisor.next()
        supervisor.next()

        slopes = supervisor.rtc.get_slopes(0)

        modes_DM0 = np.dot(S2M_DM0,slopes)
        modes_DM1 = np.dot(S2M_DM1,slopes)

        res_DM0_0[i] = modes_DM0[0]
        res_DM1_0[i] = -modes_DM1[0]

        voltage_DM1 = -M2V_DM1[:,0]*amp[i]
        voltage = np.concatenate((np.zeros(M2V_DM0.shape[0]), voltage_DM1), axis=0)
        supervisor.rtc.set_command(0, voltage)

        supervisor.next()
        supervisor.next()

        slopes = supervisor.rtc.get_slopes(0)

        modes_DM0 = np.dot(S2M_DM0,slopes)
        modes_DM1 = np.dot(S2M_DM1,slopes)

        res_DM0_1[i] = modes_DM0[0]
        res_DM1_1[i] = -modes_DM1[0]



               


    plt.figure()
    plt.plot(amp,res_DM0_0,label='LODM')
    plt.plot(amp,res_DM1_0,label='HODM')
    plt.legend()

    plt.figure()
    plt.plot(amp,res_DM0_1,label='LODM')
    plt.plot(amp,res_DM1_1,label='HODM')
    plt.legend()

    plt.show()


    pfits.writeto("../data2/res_DM1_integrating_DM0.fits", res_DM0, overwrite = True)
    # pfits.writeto("../data2/res_DM1.fits", res_DM1, overwrite = True)

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename
        embed(basename(__file__), locals())