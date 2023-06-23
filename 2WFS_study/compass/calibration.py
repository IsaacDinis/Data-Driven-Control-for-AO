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
    
    lambd = supervisor.config.p_targets[0].get_Lambda()*1e-6
    D = supervisor.config.p_tel.get_diam()
    RASC = 206265
 
    # ampli = 50
    ampli = lambd/D/20*RASC
    slopes = supervisor.rtc.get_slopes(0)
    M2S = np.zeros((slopes.shape[0], 2))

    supervisor.atmos.enable_atmos(False)

    #-----------------------------------------------
    # compute the command matrix [nmodes , nslopes]
    #-----------------------------------------------
    command = np.array([ampli,0])
    supervisor.rtc.set_command(0, command) 
    supervisor.next()
    supervisor.next()
    slopes = supervisor.rtc.get_slopes(0)/ampli
    M2S[:,0] = slopes.copy()

    command = np.array([0, ampli])
    supervisor.rtc.set_command(0, command) 
    supervisor.next()
    supervisor.next()
    slopes = supervisor.rtc.get_slopes(0)/ampli
    M2S[:,1] = slopes.copy()

    S2M = np.linalg.pinv(M2S) # [nmodes , nslopes]


    pfits.writeto('calib_mat/S2M.fits', S2M, overwrite = True)
    pfits.writeto('calib_mat/M2S.fits', M2S, overwrite = True)

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
