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
import DM_dyn
import controller

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

    bool_dm_dyn = True



    supervisor = Supervisor(config)
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(False)

    pupil = supervisor.get_s_pupil()

    n_act_DM0 = supervisor.config.p_dms[0].get_ntotact()
    n_act_DM1 = supervisor.config.p_dms[1].get_ntotact()
    n_act_bump = supervisor.config.p_dms[2].get_ntotact()


    S2M_DM0 = pfits.getdata('calib_mat/S2M_DM0.fits')
    S2M_DM1 = pfits.getdata('calib_mat/S2M_DM1.fits')
    M2V_DM0 = pfits.getdata('calib_mat/M2V_DM0.fits')
    M2V_DM1 = pfits.getdata('calib_mat/M2V_DM1.fits')

    a = np.array([1,-1]) 
    b = np.array([0.3,0])
    DM1_K = controller.K(1,a,b,S2M_DM1,M2V_DM1,stroke = np.inf)

    if bool_dm_dyn:
        HODM_dyn = DM_dyn.DM_dyn(n_act_DM1)

    voltage_DM1 = M2V_DM1[:,0] 
    voltage_DM0 = np.zeros(n_act_DM0)
    voltage_bump = np.zeros(n_act_bump)


    HODM_dead_act = 305
    HODM_dead_act_2 = 1000 #302
    HODM_dead_act_3 = 500


    


    n_iter = 5
    m = np.zeros(n_iter)

    mode = 0

    for i in range(n_iter):
        slopes = supervisor.rtc.get_slopes(0)
        # if bool_dm_dyn:     
            # voltage_DM1 = DM1_K.update_command(slopes)+M2V_DM1[:,mode]*1
            # voltage_DM1 = HODM_dyn.update_command(M2V_DM1[:,0])
        voltage = np.concatenate((voltage_DM0, voltage_DM1,voltage_bump), axis=0)
        supervisor.rtc.set_command(0,voltage)  

        m[i] = (S2M_DM1@slopes)[mode]
        supervisor.next()

 
    

    DM1_phase = supervisor.dms.get_dm_shape(1)

    plt.figure("HODM phase")
    plt.imshow(DM1_phase)
    plt.title("HODM phase")
    cbar = plt.colorbar()
    cbar.set_label(label="[um]", size=12)

    plt.figure()
    # plt.plot(command[:n_modes_applied])
    plt.plot(m)

    plt.xlabel("iter")
    plt.ylabel("amplitude [nm]")

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
