#!/usr/bin/env python

## @package   shesha.scripts.closed_loop_tilt
## @brief     script to control one mode (tilt)
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
script to control one mode (tilt)

Usage:
  closed_loop_tilt.py <parameters_filename> [options]

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
from scipy.io import savemat, loadmat

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

    if arguments["--modes"]:
        n_modes = (int(arguments["--modes"]))
    else:
        n_modes = 1
    n_modes_dd = 1
    supervisor = Supervisor(config)
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(True) 



    a = np.array([1.,-1]) 
    b = np.array([0.5,0])

    # K_dd = loadmat('Kdd.mat')['Kdd_matrix']
    # K_dd = K_dd.reshape((int(K_dd.shape[0]/2),2,n_modes_dd),order='F')

    dist = loadmat('data/single_mode_dist_phase.mat')["data"]

    K_dd = np.zeros((1,1))
    # Load command and influence matrix
    S2M = np.load('S2M.npy')
    M2V = np.load('M2V.npy')
    P2M = np.load('P2M.npy')

    bool_int = True
    bool_dist = True
    bool_lol = False
    #------------------------------------
    # control tilt mode
    #------------------------------------
    saved_mode = 0

    res_array = np.empty((n_iter,S2M.shape[0]))
    single_mode_res = np.empty(n_iter)
    single_mode_res_phase = np.empty(n_iter)
    u = np.zeros(n_iter)

    state_mat_int = np.zeros((2,2,n_modes))
    state_mat_dd = np.zeros((K_dd.shape[0],2,n_modes_dd))

    phase_count = 0

    p_diam = config.p_geom.get_pupdiam()
    n_pix = p_diam**2

    for i in range(n_iter):

        voltage = np.zeros(M2V.shape[0])  

        slopes = supervisor.rtc.get_slopes(0)
        phase = supervisor.target.get_tar_phase(0).reshape(n_pix)

        modes = np.dot(S2M,slopes)
        modes_phase = np.dot(P2M,phase)

        state_mat_int[1:,:,:] = state_mat_int[0:-1,:,:]
        state_mat_int[0,0,:] = modes[0:n_modes]
        state_mat_int[0,1,:] = 0
        command_int = np.dot(b,state_mat_int[:,0,:]) - np.dot(a,state_mat_int[:,1,:])
        # command_int -= np.mean(command_int)  
        state_mat_int[0,1,:] = command_int

        if not bool_int:
            state_mat_dd[1:,:,:] = state_mat_dd[0:-1,:,:]
            state_mat_dd[0,0,:] = modes[:n_modes_dd]
            state_mat_dd[0,1,:] = 0
            command_dd = np.sum(np.multiply(K_dd[:,0,:],state_mat_dd[:,0,:]) - np.multiply(K_dd[:,1,:],state_mat_dd[:,1,:]),0)
            state_mat_dd[0,1,:] = command_dd
            command_int[0] = command_dd[0]
        # if i > 2:
        # command_int[0] = dist[0,i+2]
        u[i] = command_int[0] 
        voltage = -M2V[:,0:n_modes] @ command_int
        if not bool_dist:
            supervisor.rtc.set_perturbation_voltage(0, "", voltage)

        single_mode_res[i] = modes[0]
        single_mode_res_phase[i] = modes_phase[0]
        # print('u = {:.3f} e = {:.3f} \n'.format(command_int[0], modes_phase[0]))
        if i%100==0:
            strehl = supervisor.target.get_strehl(0)
            print('s.e = {:.3f} l.e = {:.3f} \n'.format(strehl[0], strehl[1]))
            
        supervisor.next()
    if bool_dist:
        savemat('data/single_mode_dist_10.mat',{"data": single_mode_res[0:]})
        # savemat('data/single_mode_dist_phase.mat',{"data": single_mode_res_phase[0:]})

    elif bool_int:
        savemat('data/single_mode_res_int.mat',{"data": single_mode_res[1:]})
        savemat('data/single_mode_command_int.mat',{"data": u[1:]})
    elif bool_lol:
        savemat('data/single_mode_res_lol.mat',{"data": single_mode_res[1:]})
        savemat('data/single_mode_command_lol.mat',{"data": u[1:]})
    else:
        savemat('data/single_mode_res_dd.mat',{"data": single_mode_res[1:]})
        savemat('data/single_mode_command_dd.mat',{"data": u[1:]})

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename
        embed(basename(__file__), locals())