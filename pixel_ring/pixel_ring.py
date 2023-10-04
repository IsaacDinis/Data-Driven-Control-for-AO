# supervisor.config.p_dms[1].get_ntotact()
# supervisor.rtc.set_command(0,a)
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
from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 

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

    n_act = supervisor.config.p_dms[0].get_ntotact()
    command = np.zeros(n_act+2)


    ########################## DM tilt ###############################################
    DM_command = supervisor.config.p_dms[0].get_xpos()
    DM_command -= np.mean(DM_command)
    DM_command = DM_command/np.max(DM_command)/1.2639729 # normalisation for rms = 1 over pupil

    ########################## TT mirror tilt #########################################
    TT_command = 1/9.804306 # normalisation for rms = 1 over pupil

    ########################## PHASE tilt #########################################
    m_pupil_size = supervisor.get_m_pupil().shape[0]
    pupil_grid = make_pupil_grid(m_pupil_size)
    zernike = make_zernike_basis(3, 1, pupil_grid)
    tilt_phase = zernike[2].shaped
    tilt_record = np.dstack([tilt_phase])

    ########################## PLOT 1 #########################################

    command[:n_act] = DM_command
    command[n_act] = -TT_command
    supervisor.rtc.set_command(0,command)
    supervisor.next()
    supervisor.next()
    phase = supervisor.target.get_tar_phase(0,pupil = True)
    plt.subplots()
    plt.imshow(phase)
    plt.colorbar()
    plt.title("diff DM TT mirror")

    ########################## PLOT 2 #########################################
    supervisor.tel.set_input_phase(tilt_record)
    command[:n_act] = -DM_command
    command[n_act] = 0
    supervisor.rtc.set_command(0,command)
    supervisor.next()
    supervisor.next()
    phase = supervisor.target.get_tar_phase(0,pupil = True)
    plt.subplots()
    plt.imshow(phase)
    plt.colorbar()
    plt.title("diff DM phase tilt")


    ########################## PLOT 3 #########################################
    supervisor.tel.set_input_phase(tilt_record)
    command[:n_act] = 0
    command[n_act] = -TT_command
    supervisor.rtc.set_command(0,command)
    supervisor.next()
    supervisor.next()
    phase = supervisor.target.get_tar_phase(0,pupil = True)
    plt.subplots()
    plt.imshow(phase)
    plt.colorbar()
    plt.title("diff TT mirror phase tilt")



    plt.show()

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())