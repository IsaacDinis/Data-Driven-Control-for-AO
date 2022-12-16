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
    
    inf_mat, _ = supervisor.basis.compute_modes_to_volts_basis("KL2V") # or "KL2V" [nvolts ,nmodes]
    supervisor.rtc.set_perturbation_voltage(0, "", inf_mat[:,0])

    nactu = inf_mat.shape[0]
    nmodes = inf_mat.shape[1]
    # nmodes = 100
    ampli = 0.1
    # ampli = 0.01
    slopes = supervisor.rtc.get_slopes(0)
    imat = np.zeros((slopes.shape[0], nmodes))
    # imat = np.zeros((slopes.shape[0], nmodes+2))
    supervisor.atmos.enable_atmos(False)

    #-----------------------------------------------
    # compute the command matrix [nmodes , nslopes]
    #-----------------------------------------------
    for mode in range(nmodes):
        supervisor.rtc.set_perturbation_voltage(0, "", inf_mat[:,mode]*ampli) 
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/ampli
        imat[:,mode] = slopes.copy()

    # #tip
    # supervisor.rtc.set_perturbation_voltage(0, "", inf_mat[:,-2]*ampli) 
    # supervisor.next()
    # supervisor.next()
    # slopes = supervisor.rtc.get_slopes(0)/ampli
    # imat[:,-2] = slopes.copy()

    # #tilt
    # supervisor.rtc.set_perturbation_voltage(0, "", inf_mat[:,-1]*ampli) 
    # supervisor.next()
    # supervisor.next()
    # slopes = supervisor.rtc.get_slopes(0)/ampli
    # imat[:,-1] = slopes.copy()

    command_mat = np.linalg.pinv(imat) # [nmodes , nslopes]

    # np.savetxt('command_mat_KL2V.csv', command_mat, delimiter=",")
    # np.savetxt('inf_mat_KL2V.csv', inf_mat, delimiter=",")

    np.save('command_mat_KL2V.npy', command_mat)
    np.save('inf_mat_KL2V.npy', inf_mat)

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
