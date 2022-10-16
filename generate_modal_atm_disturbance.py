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
        n_iter = 10000

    if arguments["--modes"]:
        n_modes = (int(arguments["--modes"]))
    else:
        n_modes = 1200

    supervisor = Supervisor(config)
    supervisor.rtc.open_loop(0) # disable implemented controller

    supervisor.atmos.enable_atmos(True)

    command_mat = np.genfromtxt('../../control_matrices/command_mat_KL2V.csv',delimiter=",")
    inf_mat = np.genfromtxt('../../control_matrices/inf_mat_KL2V.csv',delimiter=",")
    
    state_mat_int = np.zeros((2,2,n_modes))
    state_mat_tt_int = np.zeros((2,2,2))

    a = np.array([1.,-1]) 
    b = np.array([0.50,0])

    command_prev = np.zeros(n_modes)
    command_prev_tt = np.zeros(2)
    #------------------------------------
    # apply disturbance
    #------------------------------------
    dist_matrix = np.empty((n_iter,n_modes+2))

    for i in range(n_iter):
        slopes = supervisor.rtc.get_slopes(0) # prendre slopes
        modes = np.dot(command_mat,slopes)

        state_mat_int[1:,:,:] = state_mat_int[0:-1,:,:]
        state_mat_int[0,0,:] = modes[0:n_modes]
        state_mat_int[0,1,:] = 0
        command_int = np.dot(b,state_mat_int[:,0,:]) - np.dot(a,state_mat_int[:,1,:])
        state_mat_int[0,1,:] = command_int 

        voltage = -inf_mat[:,0:n_modes] @ command_int

        state_mat_tt_int[1:,:,:] = state_mat_tt_int[0:-1,:,:]
        state_mat_tt_int[0,0,:] = modes[-2:]
        state_mat_tt_int[0,1,:] = 0
        command_int_tt = np.dot(b,state_mat_tt_int[:,0,:]) - np.dot(a,state_mat_tt_int[:,1,:])
        state_mat_tt_int[0,1,:] = command_int_tt

        
        voltage -= inf_mat[:,-2:] @ command_int_tt

        # dist_matrix[i,:] = np.block([modes[0:n_modes]+command_prev,modes[-2:]+command_prev_tt])
        dist_matrix[i,:] = np.block([modes[0:n_modes],modes[-2:]])
        command_prev_tt = state_mat_tt_int[1,1,:] 
        command_prev = state_mat_int[1,1,:] 

        # supervisor.rtc.set_perturbation_voltage(0, "", voltage)

        if i%100==0:
            strehl = supervisor.target.get_strehl(0)
            print('s.e = {:.5f} l.e = {:.5f} \n'.format(strehl[0], strehl[1]))

        supervisor.next()


    np.savetxt('../../disturbances/dist_matrix_'+ str(n_modes) +'_modes_KL2V_34.csv', dist_matrix, delimiter=",")
    
    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename
        embed(basename(__file__), locals())


