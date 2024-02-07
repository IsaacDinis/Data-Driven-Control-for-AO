"""
script to control one mode (tilt)

Usage:
  closed_loop_tilt.py <parameters_filename> [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
  -i, --interactive  keep the script interactive
  -d, --devices devices      Specify the devices
"""

from shesha.config import ParamConfig
from docopt import docopt
import numpy as np
from scipy.io import savemat, loadmat
import astropy.io.fits as pfits
from matplotlib import pyplot as plt
import os
from datetime import datetime
import utils
from scipy.spatial import KDTree

if __name__ == "__main__":
    arguments = docopt(__doc__)

    param_file = arguments["<parameters_filename>"]

    config = ParamConfig(param_file)

    from shesha.supervisor.compassSupervisor import CompassSupervisor as Supervisor

    if arguments["--devices"]:
        config.p_loop.set_devices([
                int(device) for device in arguments["--devices"].split(",")
        ])
    supervisor = Supervisor(config)
    
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(False) 
    p_geom = supervisor.config.p_geom

    p_dm = supervisor.config.p_dms[1]
    xpos0 = p_dm._xpos-p_geom.get_p1()
    ypos0 = p_dm._ypos-p_geom.get_p1()

    # p_dm1 = supervisor.config.p_dms[2]
    # xpos1 = p_dm1._xpos-p_geom.get_p1()
    # ypos1 = p_dm1._ypos-p_geom.get_p1()

    # act_pos = pfits.getdata("act_pos.fits")
    # act_pos *= (np.max(ypos0)-np.min(ypos0))/(np.max(act_pos[:,0])-np.min(act_pos[:,0]))
    # act_pos[:,0] += np.min(xpos0)-np.min(act_pos[:,0])
    # act_pos[:,1] += np.min(ypos0)-np.min(act_pos[:,1])


    # plt.figure()
    # plt.scatter(act_pos[:,0], act_pos[:,1], marker='.', color="blue")
    # plt.scatter(xpos0, ypos0, marker='.', color="red")
    # plt.scatter(xpos1, ypos1, marker='.', color="green")
    # plt.scatter(act_pos[383,0], act_pos[383,1], marker='.', color="yellow")
    # n = 305
    # plt.scatter(xpos1[n], ypos1[n], marker='.', color="cyan")

    n_actus_DM0 = supervisor.config.p_dms[0].get_ntotact()
    n_actus_DM1 = supervisor.config.p_dms[1].get_ntotact()

    HODM_act = 305
    pos_LODM = np.array([supervisor.config.p_dms[0].get_xpos(),supervisor.config.p_dms[0].get_ypos()]).T
    pos_HODM = np.array([supervisor.config.p_dms[1].get_xpos(),supervisor.config.p_dms[1].get_ypos()]).T
    kd_tree_LODM = KDTree(pos_LODM)
    V_DM0_2_V_DM1 = pfits.getdata('calib_mat/V_DM0_2_V_DM1.fits')

    d, i = kd_tree_LODM.query(pos_HODM[HODM_act,:], k=4)
    w = 1/d
    w /= np.sum(w)

    command_LODM = np.zeros(n_actus_DM0)
    command_dead_act = np.zeros(n_actus_DM0 + n_actus_DM1)
    for act in range(4):
        command_LODM[i[act]] = w[act]/0.1
        command_HODM = -V_DM0_2_V_DM1@command_LODM
        command_HODM[command_HODM>-0.0001] = 0
        command_dead_act += np.concatenate([command_LODM,command_HODM])
        command_LODM *= 0
    command_dead_act *= -1
    command_dead_act[n_actus_DM0+HODM_act] = 0
    # command_dead_act*= 0
    command_dead_act[n_actus_DM0+HODM_act] = 3.5
    supervisor.rtc.set_command(0, command_dead_act)
    supervisor.next()
    supervisor.next()
    supervisor.next()

    DM0_phase = supervisor.dms.get_dm_shape(0)
    DM1_phase = supervisor.dms.get_dm_shape(1)#*pupil_valid_DM1
    # DM2_phase = supervisor.dms.get_dm_shape(2)#*pupil_valid_DM1
    target_phase = supervisor.target.get_tar_phase(0,pupil=True)

    plt.figure(1)
    plt.imshow(DM0_phase)
    plt.figure(2)
    plt.imshow(DM1_phase)
    plt.figure(3)
    plt.imshow(target_phase)
    # pupil = supervisor.config.p_geom.get_spupil()
    # plt.figure()
    # plt.imshow(pupil)
    

    # np.max(act_pos[:,1])-np.min(act_pos[:,1])
    # np.max(act_pos[:,0])-np.min(act_pos[:,0])
    # np.max(ypos0)-np.min(ypos0)
    # np.max(xpos0)-np.min(xpos0)

