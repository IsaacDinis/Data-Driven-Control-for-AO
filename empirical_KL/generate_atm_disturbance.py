#!/usr/bin/env python

## @package   shesha.scripts.generate_modal_atm_disturbance
## @brief     script to generate disturbance
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
script to generate disturbance

Usage:
  generate_modal_atm_disturbance.py <parameters_filename> [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
  -i, --interactive  keep the script interactive
  -d, --devices devices      Specify the devices
  -n, --niter niter       Number of iterations
  -m, --modes modes       Number of modes to save
"""
from shesha.config import ParamConfig
from docopt import docopt
import numpy as np

if __name__ == "__main__":
    arguments = docopt(__doc__)

    param_file = arguments["<parameters_filename>"]

    config = ParamConfig(param_file)

    from shesha.supervisor.compassSupervisor import CompassSupervisor as Supervisor

    if arguments["--devices"]:
        config.p_loop.set_devices([
                int(device) for device in arguments["--devices"].split(",")
        ])

    if arguments["--niter"]:
        n_iter = (int(arguments["--niter"]))
    else:
        n_iter = 1000

    supervisor = Supervisor(config)

    n_act_eff = config.p_dms[0].get_xpos().shape[0]
    pupil_size = config.p_geom.pupdiam

    dist_matrix = np.zeros((pupil_size**2,n_iter))
    dist_matrix_act = np.zeros((n_act_eff,n_iter))

    inf_mat = np.load('inf_mat.npy')
    P2A = np.linalg.pinv(inf_mat)



    supervisor.rtc.open_loop(0) # disable implemented controller

    supervisor.atmos.enable_atmos(True)


    #------------------------------------
    # apply disturbance
    #------------------------------------
    for i in range(n_iter):
        # config.p_atmos.set_winddir([np.random.randint(0,359)]) 
        # for j in range(200):
        supervisor.next()
        phase = supervisor.target.get_tar_phase(0).reshape(pupil_size**2)
        phase -= np.mean(phase)
        dist_matrix[:,i] = phase
        dist_matrix_act[:,i] = P2A @ phase
        

    # np.save('dist_matrix.npy',dist_matrix)
    np.save('dist_matrix_act.npy',dist_matrix_act)

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename
        embed(basename(__file__), locals())


