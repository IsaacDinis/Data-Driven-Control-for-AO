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
        n_iter = 1000*1

    if arguments["--modes"]:
        n_modes = (int(arguments["--modes"]))
    else:
        n_modes = 880
    n_modes_dd = 880
    supervisor = Supervisor(config)
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(True) 

    # Load controller coefficients
    # K_tt = np.genfromtxt('../../controllers/Kdd_tt_KL2V_8.csv',delimiter=',')
    # K_DM = np.genfromtxt('../../controllers/Kdd_1st_DM_KL2V_8.csv',delimiter=',')
    K_tt = np.genfromtxt('Kdd_tt_KL2V_8.csv',delimiter=',')
    K_DM = np.genfromtxt('Kdd_1st_DM_KL2V_8.csv',delimiter=',')

    # K = K.reshape((int(K.shape[0]/2),2,1307),order='F')
    K_tt = K_tt.reshape((int(K_tt.shape[0]/2),2,2),order='F')
    K_DM = K_DM.reshape((int(K_DM.shape[0]/2),2,n_modes_dd),order='F')
    a = np.array([1.,-1]) 
    b = np.array([0.20,0])
    b_tt = np.array([0.50,0])
    # Load command and influence matrix
    # command_mat = np.genfromtxt('../../control_matrices/command_mat_KL2V.csv',delimiter=",")
    # inf_mat = np.genfromtxt('../../control_matrices/inf_mat_KL2V.csv',delimiter=",")
    command_mat = np.genfromtxt('command_mat_KL2V.csv',delimiter=",")
    inf_mat = np.genfromtxt('inf_mat_KL2V.csv',delimiter=",")

    bool_int = False
    bool_int_tt = True
    #------------------------------------
    # control tilt mode
    #------------------------------------
    saved_mode = 0
    # res_array = np.empty((n_iter,n_modes+2))
    res_array = np.empty((n_iter,command_mat.shape[0]))
    u = np.zeros(n_iter)
    strehl_se = np.zeros(n_iter)
    slopes_rms = np.zeros(n_iter)
    strehl_le = np.zeros(n_iter)
    strehl_phase = np.zeros(int(n_iter/100))

    state_mat_int = np.zeros((2,2,n_modes))
    state_mat_dd = np.zeros((K_DM.shape[0],2,n_modes_dd))
    # state_mat = np.zeros((2,2,n_modes))

    state_mat_tt_int = np.zeros((2,2,2))
    state_mat_tt_dd = np.zeros((K_tt.shape[0],2,2))
    phase_count = 0
    for i in range(n_iter):
        # if i == 1500:
        #     supervisor.target.reset_strehl(0)
        voltage = np.zeros(inf_mat.shape[0])  

        slopes = supervisor.rtc.get_slopes(0)

        modes = np.dot(command_mat,slopes)
        # state_mat[1:,:,0] = state_mat[0:-1,:,0]
        # state_mat[0,0,0] = modes[0+1]
        # state_mat[0,1,0] = 0
        # command = np.dot(K[:,0,0],state_mat[:,0,0]) - np.dot(K[:,1,0],state_mat[:,1,0])
        # state_mat[0,1,0] = command
        
        
        state_mat_int[1:,:,:] = state_mat_int[0:-1,:,:]
        state_mat_int[0,0,:] = modes[0:n_modes]
        state_mat_int[0,1,:] = 0
        command_int = np.dot(b,state_mat_int[:,0,:]) - np.dot(a,state_mat_int[:,1,:])
        state_mat_int[0,1,:] = command_int

        if not bool_int:
            state_mat_dd[1:,:,:] = state_mat_dd[0:-1,:,:]
            state_mat_dd[0,0,:] = modes[:n_modes_dd]
            state_mat_dd[0,1,:] = 0
            command_dd = np.sum(np.multiply(K_DM[:,0,:],state_mat_dd[:,0,:]) - np.multiply(K_DM[:,1,:],state_mat_dd[:,1,:]),0)
            state_mat_dd[0,1,:] = command_dd
            command_int[0:880] = command_dd[0:880]
            # print(command_dd[0])
            # res_array[i,mode] = modes[mode]
        voltage = -inf_mat[:,0:n_modes] @ command_int

        state_mat_tt_int[1:,:,:] = state_mat_tt_int[0:-1,:,:]
        state_mat_tt_int[0,0,:] = modes[-2:]
        state_mat_tt_int[0,1,:] = 0
        command_int_tt = np.dot(b,state_mat_tt_int[:,0,:]) - np.dot(a,state_mat_tt_int[:,1,:])
        state_mat_tt_int[0,1,:] = command_int_tt

        if not bool_int_tt:
            state_mat_tt_dd[1:,:,:] = state_mat_tt_dd[0:-1,:,:]
            state_mat_tt_dd[0,0,:] = modes[-2:]
            state_mat_tt_dd[0,1,:] = 0
            command_tt_dd = np.sum(np.multiply(K_tt[:,0,:],state_mat_tt_dd[:,0,:]) - np.multiply(K_tt[:,1,:],state_mat_tt_dd[:,1,:]),0)
            state_mat_tt_dd[0,1,:] = command_tt_dd
            command_int_tt = command_tt_dd

        voltage -= inf_mat[:,-2:] @ command_int_tt

        # print('input = {:.5f} command = {:.5f} \n'.format(modes[1], command))
        # if np.amax(np.absolute(voltage)) > 2:
        #     print("Warning saturation")
        #     print(np.amax(np.absolute(voltage)))
        # supervisor.rtc.set_perturbation_voltage(0, "", voltage)
        res_array[i,:] = modes
        # supervisor.target.comp_strehl(0)
        
        # res_array[i] = np.sum(modes[2:10]**2)
        # res_array[i] = np.sum(modes**2)
        strehl = supervisor.target.get_strehl(0)
        strehl_se[i] = strehl[0];
        strehl_le[i] = strehl[1];
        

        slopes_rms[i] = np.std(slopes); 

        if i%100==0:
            print('s.e = {:.5f} l.e = {:.5f} \n'.format(strehl[0], strehl[1]))
            wfs_phase = supervisor.wfs.get_wfs_phase(0)
            tar_phase = supervisor.target.get_tar_phase(0)
            # np.savetxt("phase_dd/phase_dd_"+str(phase_count)+".csv", wfs_phase, delimiter=",")
            # np.savetxt("phase_dd/tar_dd_"+str(phase_count)+".csv", tar_phase, delimiter=",")
            strehl_phase[phase_count] = strehl_se[i]
            phase_count += 1
        supervisor.next()
    # np.savetxt("phase_diff_limit/phase_diff_limit_wfs.csv", supervisor.wfs.get_wfs_phase(0), delimiter=",")
    np.savetxt("phase_diff_limit/phase_turb_tar.csv", supervisor.target.get_tar_phase(0), delimiter=",")
    # np.savetxt("phase_dd/strehl_dd.csv", strehl_phase, delimiter=",")

    # psf = supervisor.target.get_tar_image(0,expo_type = "le")
    
    # wfs_phase = supervisor.wfs.get_wfs_phase(0)

    # if bool_int:
    #     np.savetxt("../../residuals/wind34_int/strehl_int.csv", np.column_stack((u,strehl_se,strehl_le)), delimiter=",")
    #     np.savetxt("../../residuals/wind34_int/residual_int.csv", res_array, delimiter=",")
    #     np.savetxt("../../residuals/wind34_int/psf_int.csv", psf, delimiter=",")
    #     np.savetxt("../../residuals/wind34_int/wfs_phase_int.csv", wfs_phase, delimiter=",")
    #     np.savetxt("../../residuals/wind34_int/slopes_rms_int.csv", slopes_rms, delimiter=",")
    # else:
    #     np.savetxt("../../residuals/wind34_int/strehl.csv", np.column_stack((u,strehl_se,strehl_le)), delimiter=",")
    #     np.savetxt("../../residuals/wind34_int/residual.csv", res_array, delimiter=",")
    #     np.savetxt("../../residuals/wind34_int/psf.csv", psf, delimiter=",")
    #     np.savetxt("../../residuals/wind34_int/wfs_phase.csv", wfs_phase, delimiter=",")
    #     np.savetxt("../../residuals/wind34_int/slopes_rms.csv", slopes_rms, delimiter=",")
    # np.savetxt("../../residuals/psf/psf_ol_r0_08.csv", psf, delimiter=",")
    # np.savetxt("../../residuals/psf/wfs_image_int_8.csv", wfs_image, delimiter=",")
    # np.savetxt("../../residuals/psf/wfs_phase_int_8.csv", wfs_phase, delimiter=",")
    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename
        embed(basename(__file__), locals())