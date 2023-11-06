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


    M2V, _ = supervisor.basis.compute_modes_to_volts_basis("KL2V") # or "KL2V" [nvolts ,nmodes]


    n_actus_DM0 = supervisor.config.p_dms[0].get_ntotact()


    n_modes_DM0 = 1000

 
    pupil_diam = supervisor.config.p_geom.get_pupdiam()
    M2P_DM0 = np.zeros((pupil_diam*pupil_diam, n_modes_DM0))

    phase_tt_DM = np.zeros((pupil_diam, pupil_diam, 4))
    ampli = 0.01
    slopes = supervisor.rtc.get_slopes(0)
    M2S_DM0 = np.zeros((slopes.shape[0], n_modes_DM0))


    # M2S = np.zeros((slopes.shape[0], nmodes+2))
    supervisor.atmos.enable_atmos(False)

    print('N act LODM = {:d} \n'.format(n_actus_DM0))

    print('N slopes = {:d} \n'.format(slopes.shape[0]))
    
    #-----------------------------------------------
    # compute the command matrix [nmodes , nslopes]
    #-----------------------------------------------
    for mode in range(n_modes_DM0):
        supervisor.rtc.set_command(0, M2V[:,mode]*ampli) 
        supervisor.next()
        supervisor.next()
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/ampli
        target_phase = np.ndarray.flatten(supervisor.target.get_tar_phase(0,pupil=True)/ampli)

        M2S_DM0[:,mode] = slopes.copy()
        M2P_DM0[:,mode] = target_phase.copy()



    S2M_DM0 = np.linalg.pinv(M2S_DM0) # [nmodes , nslopes]




    M2V_DM0 = M2V[:n_actus_DM0,:n_modes_DM0]


    V2M_DM0 = np.linalg.pinv(M2V_DM0)




    # pfits.writeto('calib_mat/P2M_DM0.fits', P2M_DM0, overwrite = True)
    # pfits.writeto('calib_mat/P2M_DM1.fits', P2M_DM1, overwrite = True)

    p_geom = supervisor.config.p_geom
    pos_HODM = np.array([supervisor.config.p_dms[0].get_xpos(),supervisor.config.p_dms[0].get_ypos()]).T
    plt.imshow(M2P_DM0[:,0].reshape((480,480)))
    plt.scatter(pos_HODM[:,0]-p_geom.get_p1(),pos_HODM[:,1]-p_geom.get_p1(), marker='.', color="blue")
    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())