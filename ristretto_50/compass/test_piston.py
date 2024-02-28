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
from scipy.spatial import KDTree
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

    bool_flat = False
    bool_hump = False
    bool_dead_act = True
    bool_dead_act_compensation = True

    supervisor = Supervisor(config)
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(False)

    pupil = supervisor.get_s_pupil()

    n_act_DM0 = supervisor.config.p_dms[0].get_ntotact()
    n_act_DM1 = supervisor.config.p_dms[1].get_ntotact()
    n_act_bump = supervisor.config.p_dms[2].get_ntotact()

    voltage_bump = np.zeros(n_act_bump)
    # voltage_DM1 = np.ones(n_act_DM1)*0
    voltage_DM1 = np.ones(n_act_DM1)*(1.5/2) 
    voltage_DM0 = np.zeros(n_act_DM0)
    max_voltage = 3.5/2
    voltage_dead_act = max_voltage

    voltage = np.concatenate((voltage_DM0, voltage_DM1,voltage_bump), axis=0)
    if bool_dead_act :
        voltage[n_act_DM0+305] = voltage_dead_act

    voltage[n_act_DM0:n_act_DM0+n_act_DM1] = np.clip(voltage[n_act_DM0:n_act_DM0+n_act_DM1],-max_voltage,max_voltage)
    supervisor.rtc.set_command(0,voltage)
    supervisor.next()
    supervisor.next()
    supervisor.next()    

    p_geom = supervisor.config.p_geom
 

    
    n_hump = n_act_bump
    hump_offset_HODM = 9
    hump_pos = np.array([[101,60,63,63,55,57,63,51],[24,34,30,20,44,41,100,103]])+hump_offset_HODM
    dead_act_pos = np.array([108,34])

 

    DM1_phase = supervisor.dms.get_dm_shape(1)
    flat = pfits.getdata('calib_mat/flat.fits')
    if bool_flat:
        supervisor.tel.set_input_phase(flat)

    pos_LODM = np.array([supervisor.config.p_dms[0].get_xpos(),supervisor.config.p_dms[0].get_ypos()]).T
    pos_HODM = np.array([supervisor.config.p_dms[1].get_xpos(),supervisor.config.p_dms[1].get_ypos()]).T
    kd_tree_LODM = KDTree(pos_LODM)
    V_DM0_2_V_DM1 = pfits.getdata('calib_mat/V_DM0_2_V_DM1.fits')
    HODM_dead_act = 305

    d, i = kd_tree_LODM.query(pos_HODM[HODM_dead_act,:], k=4)
    w = 1/d
    w /= np.sum(w)

    command_LODM = np.zeros(n_act_DM0)
    command_dead_act = np.zeros(n_act_DM0 + n_act_DM1)
    for act in range(4):
        command_LODM[i[act]] = w[act]
        command_HODM = -V_DM0_2_V_DM1@command_LODM
        command_HODM[command_HODM>-0.0001] = 0
        command_dead_act += np.concatenate([command_LODM,command_HODM])
        command_LODM *= 0
    command_dead_act *= -1
    # command_dead_act[command_dead_act>voltage_dead_act] = voltage_dead_act
    command_dead_act *= voltage_dead_act/command_dead_act[n_act_DM0+HODM_dead_act]
    # command_dead_act *= DM1_phase[dead_act_pos[0],dead_act_pos[1]]/command_dead_act[n_act_DM0+HODM_dead_act]
    command_dead_act[n_act_DM0+HODM_dead_act] = 0
    command_dead_act = np.concatenate((command_dead_act,np.zeros(n_act_bump)), axis=0)
    
    print(DM1_phase[dead_act_pos[0],dead_act_pos[1]])

    if bool_dead_act_compensation:
        voltage += command_dead_act

   
    voltage[n_act_DM0:n_act_DM0+n_act_DM1] = np.clip(voltage[n_act_DM0:n_act_DM0+n_act_DM1],-max_voltage,max_voltage)
    supervisor.rtc.set_command(0,voltage)
    supervisor.next()
    supervisor.next()
    supervisor.next()    

    DM1_phase = supervisor.dms.get_dm_shape(1)
  
    if bool_hump:
        if DM1_phase[tuple(hump_pos[:,0])] < 0.85 :
                voltage_bump[0] =  -DM1_phase[tuple(hump_pos[:,0])]+0.85
        else:
            voltage_bump[j] = 0
        for j in range(1,n_hump):
            if  DM1_phase[tuple(hump_pos[:,j])] < -0.2 :
                voltage_bump[j] =  -DM1_phase[tuple(hump_pos[:,j])]-0.2
            else:
                voltage_bump[j] = 0
    else:
        voltage_bump *= 0

    voltage[n_act_DM0+n_act_DM1:] += voltage_bump

    voltage[n_act_DM0:n_act_DM0+n_act_DM1] = np.clip(voltage[n_act_DM0:n_act_DM0+n_act_DM1],-max_voltage,max_voltage)
    supervisor.rtc.set_command(0,voltage)
    supervisor.next()
    supervisor.next()
    supervisor.next()    

    DM0_phase = supervisor.dms.get_dm_shape(0)
    DM1_phase = supervisor.dms.get_dm_shape(1)
    DM_bump_phase = supervisor.dms.get_dm_shape(2)

    target_phase = supervisor.target.get_tar_phase(0,pupil=True)

    plt.figure()
    plt.imshow(DM0_phase)
    cbar = plt.colorbar()
    cbar.set_label(label="[um]", size=12)
    # plt.savefig(save_path+'mean_target_phase_res.png')

    plt.figure()
    plt.imshow(DM1_phase)
    cbar = plt.colorbar()
    cbar.set_label(label="[um]", size=12)
    # plt.savefig(save_path+'mean_target_phase_res.png')

    plt.figure()
    plt.imshow(DM_bump_phase)
    cbar = plt.colorbar()
    cbar.set_label(label="[um]", size=12)
    # plt.savefig(save_path+'mean_target_phase_res.png')

    plt.figure()
    plt.imshow(target_phase)
    cbar = plt.colorbar()
    cbar.set_label(label="[um]", size=12)
    # plt.savefig(save_path+'mean_target_phase_res.png')

    # print(np.unravel_index(target_phase.argmax(), target_phase.shape))
    # print(np.unravel_index(DM1_phase.argmax(), DM1_phase.shape))
    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
