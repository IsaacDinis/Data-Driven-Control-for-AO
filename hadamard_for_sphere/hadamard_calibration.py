#!/usr/bin/env python

## @brief     hadamard calibration in compass for saxo
## @author    Isaac Dinis <https://github.com/IsaacDinis>
## @date      20.12.2022
## @copyright GNU Lesser General Public License

"""
script to generate slopes to actuators amptude (volts) and vice versa matrices using hadamard calibration method. 
cf "Fast calibration of high-order adaptive optics systems", Kasper and al. 

Usage:
  hadamard_calibration.py <parameters_filename> [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h  --help          Show this help message and exit
  -i, --interactive  keep the script interactive
  -d, --devices devices      Specify the devices
"""
from shesha.config import ParamConfig
from docopt import docopt
import numpy as np
from scipy.linalg import hadamard

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
    supervisor.atmos.enable_atmos(False) # disable atmosphere turbulence
    
    n_act = config.p_dms[0].get_xpos().shape[0] # get number of visible actuators (1377 with SAXO HODM)
    n_slopes = supervisor.rtc.get_slopes(0).shape[0] # get number of slopes
 
    amp = 0.1 # calibration voltage

    hadamard_size = 2**n_act.bit_length() # find nearest power of two
    hadamard_matrix = hadamard(hadamard_size) 

    C = np.zeros((n_slopes,hadamard_size))

    #-----------------------------------------------
    # do calibration
    #-----------------------------------------------
    for i in range(hadamard_size):
        supervisor.rtc.set_perturbation_voltage(0, "", hadamard_matrix[:n_act,i]*amp)
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/amp
        C[:,i] = slopes

    D = C @ np.linalg.inv(hadamard_matrix)
    D = D[:,:n_act]
    S2A = np.linalg.pinv(D)

    #-----------------------------------------------
    # save matrices
    #-----------------------------------------------
    np.save('S2A.npy', S2A) # slopes to amplitude
    np.save('A2S.npy', D) # amplitude to slopes

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
