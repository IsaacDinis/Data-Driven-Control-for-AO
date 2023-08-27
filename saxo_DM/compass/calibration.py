# wao.supervisor.config.p_dms[1].get_ntotact()
# wao.supervisor.rtc.set_command(0,a)
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
from shesha.util.slopesCovariance import KLmodes

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
    
    # M2V, _ = supervisor.basis.compute_modes_to_volts_basis("KL2V") # or "KL2V" [nvolts ,nmodes]

    xpos0 = supervisor.config.p_dms[0]._xpos # actus positions
    ypos0 = supervisor.config.p_dms[0]._ypos

    xpos1 = supervisor.config.p_dms[1]._xpos # actus positions
    ypos1 = supervisor.config.p_dms[1]._ypos
    L0 = 25  # [m]
    B0, l = KLmodes(xpos0, ypos0, L0, True) #basis on saxo stage
    B1, l = KLmodes(xpos1, ypos1, L0, True) #basis on saxoplus stage
    M2V = np.zeros((B0.shape[0]+B1.shape[0],B0.shape[1]+B1.shape[1]))
    M2V[:B0.shape[0],:B0.shape[1]] = B0
    M2V[B0.shape[0]:,B0.shape[1]:] = B1


    n_actus_DM0 = supervisor.config.p_dms[0].get_ntotact()
    n_actus_DM1 = supervisor.config.p_dms[1].get_ntotact()

    n_modes_DM0 = 800
    # n_modes_DM1 = n_actus_DM1
    n_modes_DM1 = 200
 
    pupil_diam = supervisor.config.p_geom.get_pupdiam()
    phase_mode_DMO = np.zeros((pupil_diam, pupil_diam, n_modes_DM0))
    phase_mode_DM1 = np.zeros((pupil_diam, pupil_diam, n_modes_DM1))
    phase_tt_DM = np.zeros((pupil_diam, pupil_diam, 4))
    # ampli = 50
    ampli = 0.01
    slopes = supervisor.rtc.get_slopes(0)
    M2S_DM0 = np.zeros((slopes.shape[0], n_modes_DM0))
    M2S_DM1 = np.zeros((slopes.shape[0], n_modes_DM1))
    modes_stage1_2_modes_stage2 = np.zeros((n_modes_DM1, n_modes_DM1))
    # M2S = np.zeros((slopes.shape[0], nmodes+2))
    supervisor.atmos.enable_atmos(False)

    #-----------------------------------------------
    # compute the command matrix [nmodes , nslopes]
    #-----------------------------------------------

    for mode in range(n_modes_DM1):
        supervisor.rtc.set_command(0, M2V[:,n_actus_DM0+mode]*ampli) 
        supervisor.next()
        supervisor.next()
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/ampli
        target_phase = supervisor.target.get_tar_phase(0,pupil=True)/ampli

        M2S_DM1[:,mode] = slopes.copy()
        phase_mode_DM1[:,:,mode] = target_phase.copy()

    S2M_DM1 = np.linalg.pinv(M2S_DM1) # [nmodes , nslopes

    for mode in range(n_modes_DM0):
        supervisor.rtc.set_command(0, M2V[:,mode]*ampli) 
        supervisor.next()
        supervisor.next()
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/ampli
        target_phase = supervisor.target.get_tar_phase(0,pupil=True)/ampli

        M2S_DM0[:,mode] = slopes.copy()
        phase_mode_DMO[:,:,mode] = target_phase.copy()
        if mode < n_modes_DM1:
            modes_stage1_2_modes_stage2[:,mode] = S2M_DM1@slopes  

    S2M_DM0 = np.linalg.pinv(M2S_DM0) # [nmodes , nslopes]



    
    



    M2V_DM0 = M2V[:n_actus_DM0,:n_modes_DM0]
    M2V_DM1 = M2V[n_actus_DM0:n_actus_DM0+n_actus_DM1,n_actus_DM0:n_actus_DM0+n_modes_DM1]


    n_projected_modes = 150 # number of modes projected from saxo HODM to saxoplus HODM,
    # S2M_stage1 = np.linalg.pinv(imat1[:, :n_projected_modes]) # SAXO imat
    # V2M_stage1 = np.linalg.pinv(B0[:, :n_projected_modes]) # SAXO KL basis

    S2M_stage2 = S2M_DM1
    V2M_stage1 = np.linalg.pinv(M2V_DM0[:,:n_projected_modes])

    # modes_stage1_2_modes_stage2 = S2M_stage1@imat01[:,:n_projected_modes] # Modal projection matrix
    # volts_stage1_2_volts_stage2 = B1[:, :n_projected_modes]@modes_stage1_2_modes_stage2@V2M_stage1 # Voltage projection matrix

    modes_stage1_2_modes_stage2 = S2M_stage2[:n_projected_modes,:]@M2S_DM0[:,:n_projected_modes] # Modal projection matrix
    volts_stage1_2_volts_stage2 = M2V_DM1[:, :n_projected_modes]@modes_stage1_2_modes_stage2@V2M_stage1 # Voltage projection matrix



    ################## Used for integrated solution #####################################



    V2M_DM0 = np.linalg.pinv(M2V_DM0)

    M_DM0_2_M_DM1 = modes_stage1_2_modes_stage2
    V_DM0_2_V_DM1 = volts_stage1_2_volts_stage2
    V_DM1_2_V_DM0 = np.linalg.pinv(V_DM0_2_V_DM1)

    np.save('calib_mat/S2M_DM0.npy', S2M_DM0)

    np.save('calib_mat/S2M_DM1.npy', S2M_DM1)
    np.save('calib_mat/M2V.npy', M2V)
    np.save('calib_mat/M2V_DM0.npy', M2V_DM0)

    np.save('calib_mat/M2V_DM1.npy', M2V_DM1)
    np.save('calib_mat/M_DM0_2_M_DM1.npy', M_DM0_2_M_DM1)
    np.save('calib_mat/V_DM0_2_V_DM1.npy', V_DM0_2_V_DM1)
    np.save('calib_mat/V_DM1_2_V_DM0.npy', V_DM1_2_V_DM0)

    pfits.writeto("../data_parallel/phase_mode_DMO.fits", phase_mode_DMO, overwrite = True)
    pfits.writeto("../data_parallel/phase_mode_DM1.fits", phase_mode_DM1, overwrite = True)

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())