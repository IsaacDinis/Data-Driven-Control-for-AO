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
import utils
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
    supervisor.atmos.enable_atmos(False)

    pupil = supervisor.get_s_pupil()

    pos_bump = np.array([supervisor.config.p_dms[2].get_ypos(),supervisor.config.p_dms[2].get_xpos()]).T
    pos_dm = np.array([supervisor.config.p_dms[1].get_ypos(),supervisor.config.p_dms[1].get_xpos()]).T
    # plt.imshow(pupil)
    # plt.scatter(pos_bump[:,0]-p_geom.get_p1(),pos_bump[:,1]-p_geom.get_p1(), marker='.', color="green")

    n_actus_DM0 = supervisor.config.p_dms[0].get_ntotact()
    n_actus_DM1 = supervisor.config.p_dms[1].get_ntotact()
    n_actus_bump = supervisor.config.p_dms[2].get_ntotact()
    p_geom = supervisor.config.p_geom
    command = supervisor.rtc.get_command(0)
    # ampli1 = -1
    # ampli2 = -1
    # ampli3 = -1
    # ampli4 = -1

    # command[n_actus_DM0:n_actus_DM0+n_actus_DM1] = 0
    # command[n_actus_DM0+376]  = ampli1
    # command[n_actus_DM0+377]  = ampli2
    # command[n_actus_DM0+416]  = ampli3
    # command[n_actus_DM0+417]  = ampli4
    # mean_ampli = (0.6*ampli1+0.4*ampli2+0.7*ampli3+0.3*ampli4)/2

    # if  mean_ampli*2 > 1.8:
    #     command[-1] = -1.95*mean_ampli+1.8
    # else:
    #     command[-1] = 0
    # DM_phase = supervisor.dms.get_dm_shape(1)
    # if  np.max(DM_phase) - DM_phase[323,162] > 1.8:
    #     command[-1] = np.max(DM_phase) - DM_phase[323,162]-1.8
    # else:
    #     command[-1] = 0

    # if  DM_phase[323,162] < -0.1:
    #     command[-1] = -DM_phase[323,162]-0.1
    # else:
    #     command[-1] = 0  
    # command[-1] = 1 
    # command[-1] = 0
    # command[-1] = -DM_phase[301,140]
    command*=0
    command[n_actus_DM0+n_actus_DM1+6]  = 1
    command[n_actus_DM0+500]  = -1
    supervisor.rtc.set_perturbation_voltage(0, "tmp", command)
    supervisor.next()
    supervisor.next()


    phase = supervisor.target.get_tar_phase(0, pupil = True)
    hump = 6
    plt.figure()
    plt.imshow(phase)
    plt.scatter(pos_bump[hump,0]-p_geom.get_p1(),pos_bump[hump,1]-p_geom.get_p1(), marker='.', color="red")
    plt.scatter(pos_dm[:,0]-p_geom.get_p1(),pos_dm[:,1]-p_geom.get_p1(), marker='.', color="green")
    # plt.scatter(pos_dm[500,0]-p_geom.get_p1(),pos_dm[500,1]-p_geom.get_p1(), marker='.', color="yellow")
    # print(phase[301,140])
    print(np.unravel_index(phase.argmax(), phase.shape))
    supervisor.config.p_dms[2]._n1
    DM1_shape = supervisor.dms.get_dm_shape(1)
    DM2_shape = supervisor.dms.get_dm_shape(2)
    plt.figure()
    plt.imshow(DM1_shape[22:-22,22:-22]+DM2_shape[8:-8,8:-8]-phase)
    print(DM2_shape[287,149])
    # plt.imshow(supervisor.dms.get_dm_shape(2)[8:-8,8:-8])
    # np.unravel_index(DM_phase.argmax(), DM_phase.shape)
    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
