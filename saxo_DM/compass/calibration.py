# wao.supervisor.config.p_dms[1].get_ntotact()
# wao.supervisor.rtc.set_command(0,a)
# ipython -i shesha/widgets/widget_ao.py ~/Data-Driven-Control-for-AO/2DM_study/compass/compass_param.py 

"""
script to generate the command matrix

Usage:
  generate_command_matrix.py <parameters_filename> [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h  --help          Show this help message and exit
  -i, --interactive  keep the script interactive
  -d, --devices devices      Specify the devices
"""

from shesha.config import ParamConfig
from docopt import docopt
import numpy as np
import astropy.io.fits as pfits

if __name__ == "__main__":
    arguments = docopt(__doc__)

    param_file = arguments["<parameters_filename>"]

    # Get parameters from file
    config = ParamConfig(param_file)

    from shesha.supervisor.compassSupervisor import CompassSupervisor as Supervisor

    if arguments["--devices"]:
        config.p_loop.set_devices([
                int(device) for device in arguments["--devices"].split(",")
        ])

    supervisor = Supervisor(config)
    supervisor.rtc.open_loop(0) # disable implemented controller
    

    # norm = np.linalg.norm(M2V, axis = 0)
    # M2V /= norm



    n_actus_DM0 = 2
    n_actus_DM1 = 2

    n_modes_DM0 = 2
    # n_modes_DM1 = n_actus_DM1
    n_modes_DM1 = 2
 

 
    # ampli = 50
    ampli = 0.01
    slopes = supervisor.rtc.get_slopes(0)
    M2S_DM0 = np.zeros((slopes.shape[0], n_modes_DM0))
    M2S_DM1 = np.zeros((slopes.shape[0], n_modes_DM1))

    # M2S = np.zeros((slopes.shape[0], nmodes+2))
    supervisor.atmos.enable_atmos(False)

    #-----------------------------------------------
    # compute the command matrix [nmodes , nslopes]
    #-----------------------------------------------
    for mode in range(n_modes_DM0):
        supervisor.rtc.set_command(0, M2V[:,mode]*ampli) 
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/ampli
        M2S_DM0[:,mode] = slopes.copy()

    for mode in range(n_modes_DM1):
        supervisor.rtc.set_command(0, M2V[:,n_modes_DM0+mode]*ampli) 
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/ampli
        M2S_DM1[:,mode] = slopes.copy()


    S2M_DM0 = np.linalg.pinv(M2S_DM0) # [nmodes , nslopes]
    S2M_DM1 = np.linalg.pinv(M2S_DM1) # [nmodes , nslopes

    M2V_DM0 = M2V[:n_actus_DM0,:n_modes_DM0]
    M2V_DM1 = M2V[n_actus_DM0:n_actus_DM0+n_actus_DM1,n_actus_DM0:n_actus_DM0+n_modes_DM1]

    V2M_DM0 = np.linalg.pinv(M2V_DM0)

    M_DM0_2_M_DM1 = S2M_DM1@M2S_DM0
    V_DM0_2_V_DM1 = M2V_DM1@M_DM0_2_M_DM1@V2M_DM0
    V_DM1_2_V_DM0 = np.linalg.pinv(V_DM0_2_V_DM1)
    pfits.writeto("../data2/M2M.fits", M_DM0_2_M_DM1, overwrite = True)
    np.save('calib_mat/S2M_DM0.npy', S2M_DM0)
    pfits.writeto("../data3/S2M_DM0.fits", S2M_DM0, overwrite = True)
    np.save('calib_mat/S2M_DM1.npy', S2M_DM1)
    np.save('calib_mat/M2V.npy', M2V)
    np.save('calib_mat/M2V_DM0.npy', M2V_DM0)
    pfits.writeto("../data3/M2V_DM0.fits", M2V_DM0, overwrite = True)
    np.save('calib_mat/M2V_DM1.npy', M2V_DM1)
    np.save('calib_mat/M_DM0_2_M_DM1.npy', M_DM0_2_M_DM1)
    np.save('calib_mat/V_DM0_2_V_DM1.npy', V_DM0_2_V_DM1)
    np.save('calib_mat/V_DM1_2_V_DM0.npy', V_DM1_2_V_DM0)
    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
