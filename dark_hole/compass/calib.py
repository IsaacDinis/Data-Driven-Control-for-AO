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
from shesha.util.slopesCovariance import KLmodes
import astropy.io.fits as pfits
from matplotlib import pyplot as plt
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
    pupil = supervisor.get_s_pupil()

    xpos0 = supervisor.config.p_dms[0]._xpos # actus positions
    ypos0 = supervisor.config.p_dms[0]._ypos

    xpos1 = supervisor.config.p_dms[1]._xpos # actus positions
    ypos1 = supervisor.config.p_dms[1]._ypos

    L0 = 25  # [m]
    M2V_DM0, l = KLmodes(xpos0, ypos0, L0, True) #basis on saxo stage


    # norm = np.linalg.norm(M2V, axis = 0)
    # M2V /= norm

    n_actus_DM0 = supervisor.config.p_dms[0].get_ntotact()
    n_modes_DM0 = 400

    print(n_actus_DM1)

    M2V_DM0 = M2V_DM0[:,:n_modes_DM0]

    ampli = 0.1
    # ampli = 0.01
    slopes = supervisor.rtc.get_slopes(0)
    M2S_DM0 = np.zeros((slopes.shape[0], n_modes_DM0))


    # M2S = np.zeros((slopes.shape[0], n_modes+2))
    supervisor.atmos.enable_atmos(False)

    #-----------------------------------------------
    # compute the command matrix [n_modes , nslopes]
    #-----------------------------------------------
    for mode in range(n_modes_DM0):
        command = np.concatenate((M2V_DM0[:,mode]*ampli,np.zeros(n_actus_DM1+n_actus_bump)), axis=0)
        # command = np.concatenate((M2V_DM0[:,mode]*ampli,np.zeros(n_actus_DM1)), axis=0)
        supervisor.rtc.set_command(0, command) 
        supervisor.next()
        supervisor.next()
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/ampli
        M2S_DM0[:,mode] = slopes.copy()


    S2M_DM0 = np.linalg.pinv(M2S_DM0) # [n_modes , nslopes]
    V2M_DM0 = np.linalg.pinv(M2V_DM0)


    pfits.writeto('calib_mat/S2M_DM0.fits', S2M_DM0, overwrite = True)
    pfits.writeto('calib_mat/M2V_DM0.fits', M2V_DM0, overwrite = True)


 
    

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())