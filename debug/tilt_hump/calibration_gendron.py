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


    L0 = 25  # [m]
    M2V_DM0, l = KLmodes(xpos0, ypos0, L0, True) #basis on saxo stage


    n_actus_DM0 = supervisor.config.p_dms[0].get_ntotact()

    n_modes_DM0 = 1000



    M2V_DM0 = M2V_DM0[:,:n_modes_DM0]

    ampli = 0.01
    # ampli = 0.01
    slopes = supervisor.rtc.get_slopes(0)
    M2S_DM0 = np.zeros((slopes.shape[0], n_modes_DM0))


    # M2S = np.zeros((slopes.shape[0], n_modes+2))
    supervisor.atmos.enable_atmos(False)

    #-----------------------------------------------
    # compute the command matrix [n_modes , nslopes]
    #-----------------------------------------------
    for mode in range(n_modes_DM0):
        command = np.concatenate((M2V_DM0[:,mode]*ampli,np.zeros(2)), axis=0)
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

    # pfits.writeto('calib_mat/V_DM1_2_V_DM0.fits', V_DM1_2_V_DM0, overwrite = True)

    p_geom = supervisor.config.p_geom
    pos_HODM = np.array([supervisor.config.p_dms[0].get_xpos(),supervisor.config.p_dms[0].get_ypos()]).T

    plt.imshow(pupil)
    
    plt.scatter(pos_HODM[:,0]-p_geom.get_p1(),pos_HODM[:,1]-p_geom.get_p1(), marker='.', color="blue")
    plt.scatter(pos_HODM[376,0]-p_geom.get_p1(),pos_HODM[376,1]-p_geom.get_p1(), marker='.', color="red")


    
    plt.figure()

    command = np.zeros(n_actus_DM0+2)
    ampli = 1/65.47307333333335
    command[-1]  = ampli
    supervisor.rtc.set_perturbation_voltage(0, "tmp", command)
    supervisor.next()
    supervisor.next()
    slopes = supervisor.rtc.get_slopes(0)/ampli
    target_phase = supervisor.target.get_tar_phase(0,pupil=True)
    plt.imshow(target_phase)
    print(np.max(target_phase)-np.min(target_phase))

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
