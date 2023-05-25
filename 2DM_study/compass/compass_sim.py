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
    supervisor.atmos.enable_atmos(True) 

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
    res_DM0 = np.zeros(n_iter)
    res_DM1 = np.zeros(n_iter)

    #------------------------------------
    # control tilt mode
    #------------------------------------

    # res_array = np.empty((n_iter,S2M.shape[0]))
    # single_mode_res = np.empty(n_iter)

    state_mat_DM0 = np.zeros((2,2,n_modes_DM0))
    state_mat_DM1 = np.zeros((2,2,n_modes_DM1))

    bool_DMO = True
    rms_stroke = 0;
    for i in range(n_iter):


        slopes = supervisor.rtc.get_slopes(0)

        modes_DM0 = np.dot(S2M_DM0,slopes)
        modes_DM1 = np.dot(S2M_DM1,slopes)

        # if  i%4==0:
        # modes_DM0 = np.dot(S2M_DM0,slopes)
        state_mat_DM0[1:,:,:] = state_mat_DM0[0:-1,:,:]
        state_mat_DM0[0,0,:] = modes_DM0[0:n_modes_DM0]
        state_mat_DM0[0,1,:] = 0
        command_int_DM0 = np.dot(b,state_mat_DM0[:,0,:]) - np.dot(a,state_mat_DM0[:,1,:])
        # command_int -= np.mean(command_int)  
        state_mat_DM0[0,1,:] = command_int_DM0
        voltage_DM0 = -M2V_DM0[:,0:n_modes_DM0] @ command_int_DM0

        state_mat_DM1[1:,:,:] = state_mat_DM1[0:-1,:,:]
        state_mat_DM1[0,0,:] = modes_DM1[0:n_modes_DM1]
        state_mat_DM1[0,1,:] = 0
        command_int_DM1 = np.dot(b,state_mat_DM1[:,0,:]) - np.dot(a,state_mat_DM1[:,1,:])
        command_int_DM1 -= np.mean(command_int_DM1)  
        state_mat_DM1[0,1,:] = command_int_DM1
        voltage_DM1 = -M2V_DM1[:,0:n_modes_DM1] @ command_int_DM1
        voltage_DM1 *= 0

        if bool_DMO:
            # voltage_DM1 -= np.dot(V_DM0_2_V_DM1,voltage_DM0)
            voltage = np.concatenate((voltage_DM0, voltage_DM1), axis=0)
        else:
            voltage = np.concatenate((np.zeros(M2V_DM0.shape[0]), voltage_DM1), axis=0)

        supervisor.rtc.set_command(0, voltage)

        strehl = supervisor.target.get_strehl(0)

        rms_stroke += np.std(voltage_DM1)

        if i%100==0 and i > 200:
            print('s.e = {:.5f} l.e = {:.5f} \n'.format(strehl[0], strehl[1]))

        res_DM0[i] = modes_DM0[0]
        res_DM1[i] = modes_DM1[0]
        supervisor.next()
    rms_stroke /= n_iter
    print('rms_stroke = {:.5f} \n'.format(rms_stroke))

    pfits.writeto("../data2/res_DM0_alone.fits", res_DM0, overwrite = True)
    # pfits.writeto("../data2/res_DM1.fits", res_DM1, overwrite = True)

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename
        embed(basename(__file__), locals())