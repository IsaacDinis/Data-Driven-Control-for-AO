#!/usr/bin/env python

## @package   shesha.scripts.generate_command_matrix
## @brief     script to generate the command matrix
## @author    Isaac Dinis <https://github.com/IsaacDinis>
## @date      01.04.2022
## @copyright GNU Lesser General Public License
#
#  This file is part of COMPASS <https://anr-compass.github.io/compass/>
#
#  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
#  All rights reserved.
#  Distributed under GNU - LGPL
#
#  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
#  General Public License as published by the Free Software Foundation, either version 3 of the License,
#  or any later version.
#
#  COMPASS: End-to-end AO simulation tool using GPU acceleration
#  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
#
#  The final product includes a software package for simulating all the critical subcomponents of AO,
#  particularly in the context of the ELT and a real-time core inf_matd on several control approaches,
#  with performances consistent with its integration into an instrument. Taking advantage of the specific
#  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
#  conduct large simulation campaigns called to the ELT.
#
#  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
#  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
#  various systems configurations such as multi-conjugate AO.
#
#  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
#  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
#  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.

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
from scipy.linalg import hadamard

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
    

    # supervisor.rtc.set_perturbation_voltage(0, "", inf_mat[:,0])
    supervisor.atmos.enable_atmos(False)
    
    n_act = config.p_dms[0].get_xpos().shape[0] # get number of visible actuators
    
    n_slopes = supervisor.rtc.get_slopes(0).shape[0]
 
    ampli = 0.1

    hadamard_size = 2**n_act.bit_length()
    hadamard_matrix = hadamard(hadamard_size)

    C = np.zeros((n_slopes,hadamard_size))

    # # ampli = 0.01
    # slopes = supervisor.rtc.get_slopes(0)
    # imat = np.zeros((slopes.shape[0], nmodes))
    # # imat = np.zeros((slopes.shape[0], nmodes+2))
    


    #-----------------------------------------------
    # compute the command matrix [nmodes , nslopes]
    #-----------------------------------------------
    for i in range(hadamard_size):
        supervisor.rtc.set_perturbation_voltage(0, "", hadamard_matrix[:n_act,i]*ampli)
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/ampli
        C[:,i] = slopes

    D = C @ np.linalg.inv(hadamard_matrix)
    # D = D[:,:n_act_square]
    # D = D[:,:n_act_eff]
    S2A = np.linalg.pinv(D)

    # command_mat = np.linalg.pinv(imat) # [nmodes , nslopes]

    np.save('S2A.npy', S2A)
    np.save('A2S.npy', D)
    # np.savetxt('../../control_matrices/inf_mat_KL2V.csv', inf_mat, delimiter=",")

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
