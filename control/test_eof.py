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
import utils
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

    supervisor = Supervisor(config)
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(False)

    pupil = supervisor.get_s_pupil()

    n_act = supervisor.config.p_dms[0].get_ntotact()



    S2M = pfits.getdata('calib_mat/S2M.fits')
    M2V = pfits.getdata('calib_mat/M2V.fits')

    a = np.array([1,-0.99]) 
    b = np.array([0.7,0])


    # Load command and influence matrix


    S2M = pfits.getdata('calib_mat/S2M.fits')

    M2V = pfits.getdata('calib_mat/M2V.fits')

    K = controller.K(1,a,b,S2M,M2V,stroke = np.inf)
    K_eof = eof.eof(5,S2M, M2V,100)
    while(K_eof.is_trained == 0):
        slopes = supervisor.rtc.get_slopes(0)
        voltage = DM0_K.update_command(slopes)
        K_eof.train(slopes,voltage)
        supervisor.rtc.set_command
        supervisor.next()
    n_iter = 5
    m = np.zeros(n_iter)

    mode = 0
    supervisor.next()
    supervisor.next()
    for i in range(n_iter):
        
        supervisor.rtc.set_command(0,M2V[:,0])  
        slopes = supervisor.rtc.get_slopes(0)
        m[i] = (S2M@slopes)[mode]
        supervisor.next()

 
    

    DM1_phase = supervisor.dms.get_dm_shape(0)

    plt.figure("HODM phase")
    plt.imshow(DM1_phase)
    plt.title("HODM phase")
    cbar = plt.colorbar()
    cbar.set_label(label="[um]", size=12)

    plt.figure("WFS image")
    plt.imshow(supervisor.wfs.get_wfs_image(0))
    plt.title("WFS image")
    cbar = plt.colorbar()
    cbar.set_label(label="[um]", size=12)

    plt.figure()
    # plt.plot(command[:n_modes_applied])
    plt.plot(m)

    plt.xlabel("iter")
    plt.ylabel("amplitude [nm]")

    plt.show()
    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname