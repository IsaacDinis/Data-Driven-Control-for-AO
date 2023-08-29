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

    pupil_diam = supervisor.config.p_geom.get_pupdiam()
    phase_mode_tt = np.zeros((pupil_diam, pupil_diam, 2))


    # ampli = 50
    ampli = 0.1
    slopes = supervisor.rtc.get_slopes(0)
    M2S_DM0 = np.zeros((slopes.shape[0], 2))
    M2S_DM1 = np.zeros((slopes.shape[0], 2))

    # M2S = np.zeros((slopes.shape[0], nmodes+2))
    supervisor.atmos.enable_atmos(False)
    command = supervisor.rtc.get_command(0)
    command *= 0
    #-----------------------------------------------
    # compute the command matrix [nmodes , nslopes]
    #-----------------------------------------------
    for mode in range(2):
        command[mode] = ampli
        supervisor.rtc.set_command(0, command) 
        command *= 0
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/ampli
        M2S_DM0[:,mode] = slopes.copy()
        target_phase = supervisor.target.get_tar_phase(0,pupil=True)/ampli
        phase_mode_tt[:,:,mode] = target_phase.copy()

    for mode in range(2):
        command[mode+2] = ampli
        supervisor.rtc.set_command(0, command) 
        command *= 0
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/ampli
        M2S_DM1[:,mode] = slopes.copy()


    S2M_DM0 = np.linalg.pinv(M2S_DM0) # [nmodes , nslopes]
    S2M_DM1 = np.linalg.pinv(M2S_DM1) # [nmodes , nslopes

    M2V_DM0 = np.eye(2)
    M2V_DM1 = np.eye(2)

    V2M_DM0 = np.linalg.pinv(M2V_DM0)

    M_DM0_2_M_DM1 = S2M_DM1@M2S_DM0
    V_DM0_2_V_DM1 = M2V_DM1@M_DM0_2_M_DM1@V2M_DM0
    V_DM1_2_V_DM0 = np.linalg.pinv(V_DM0_2_V_DM1)
    pfits.writeto('calib_mat/S2M_DM0.fits', S2M_DM0, overwrite = True)
    pfits.writeto('calib_mat/S2M_DM1.fits', S2M_DM1, overwrite = True)
    pfits.writeto('calib_mat/M2V_DM0.fits', M2V_DM0, overwrite = True)
    pfits.writeto('calib_mat/M2V_DM1.fits', M2V_DM1, overwrite = True)
    pfits.writeto('calib_mat/M_DM0_2_M_DM1.fits', M_DM0_2_M_DM1, overwrite = True)
    pfits.writeto('calib_mat/V_DM0_2_V_DM1.fits', V_DM0_2_V_DM1, overwrite = True)
    pfits.writeto('calib_mat/V_DM1_2_V_DM0.fits', V_DM1_2_V_DM0, overwrite = True)
    pfits.writeto("../data_parallel/phase_mode_tt.fits", phase_mode_tt, overwrite = True)
    
    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
