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
    
    M2V, _ = supervisor.basis.compute_modes_to_volts_basis("KL2V") # or "KL2V" [nvolts ,nmodes]
    # norm = np.linalg.norm(M2V, axis = 0)
    # M2V /= norm



    n_actus_DM0 = supervisor.config.p_dms[0].get_ntotact()
    n_actus_DM1 = supervisor.config.p_dms[1].get_ntotact()

    n_modes_DM0 = n_actus_DM0
    # n_modes_DM1 = n_actus_DM1
    n_modes_DM1 = 800
 
    nmodes = M2V.shape[1]

    pupil_diam = supervisor.config.p_geom.get_pupdiam()
    pupil_grid = make_pupil_grid(pupil_diam)
    zernike_basis = make_zernike_basis(3, 1, pupil_grid)
    pupil_valid = zernike_basis[0].shaped
 
    # ampli = 50
    ampli = 0.01
    n_phase = int(np.sum(np.array(pupil_valid)))

    M2P_DM0 = np.zeros((n_phase, n_modes_DM0))
    M2P_DM1 = np.zeros((n_phase, n_modes_DM1))

    # M2S = np.zeros((phase.shape[0], nmodes+2))
    supervisor.atmos.enable_atmos(False)

    #-----------------------------------------------
    # compute the command matrix [nmodes , nphase]
    #-----------------------------------------------
    for mode in range(n_modes_DM0):
        supervisor.rtc.set_command(0, M2V[:,mode]*ampli) 
        supervisor.next()
        supervisor.next()
        phase = supervisor.target.get_tar_phase(0,pupil=True)/ampli
        phase = phase[pupil_valid == 1]
        M2P_DM0[:,mode] = phase.copy()

    for mode in range(n_modes_DM1):
        supervisor.rtc.set_command(0, M2V[:,n_modes_DM0+mode]*ampli) 
        supervisor.next()
        supervisor.next()
        phase = supervisor.target.get_tar_phase(0,pupil=True)/ampli
        phase = phase[pupil_valid == 1]
        M2P_DM1[:,mode] = phase.copy()


    P2M_DM0 = np.linalg.pinv(M2P_DM0) # [nmodes , nphase]
    P2M_DM1 = np.linalg.pinv(M2P_DM1) # [nmodes , nphase

    M2V_DM0 = M2V[:n_actus_DM0,:n_modes_DM0]
    M2V_DM1 = M2V[n_actus_DM0:n_actus_DM0+n_actus_DM1,n_actus_DM0:n_actus_DM0+n_modes_DM1]

    V2M_DM0 = np.linalg.pinv(M2V_DM0)

    M_DM0_2_M_DM1 = P2M_DM1@M2P_DM0
    V_DM0_2_V_DM1 = M2V_DM1@M_DM0_2_M_DM1@V2M_DM0
    V_DM1_2_V_DM0 = np.linalg.pinv(V_DM0_2_V_DM1)
    
    pfits.writeto("calib_mat/P2M_DM0.fits", P2M_DM0, overwrite = True)

    pfits.writeto("calib_mat/P2M_DM1.fits", P2M_DM1, overwrite = True)

    pfits.writeto("calib_mat/M2P_DM0.fits", M2P_DM0, overwrite = True)

    pfits.writeto("calib_mat/M2P_DM1.fits", M2P_DM1, overwrite = True)

    pfits.writeto("calib_mat/M2V.fits", M2V, overwrite = True)
 
    pfits.writeto("calib_mat/M2V_DM0.fits", M2V_DM0, overwrite = True)
    
    pfits.writeto("calib_mat/M2V_DM1.fits", M2V_DM1, overwrite = True)
    
    pfits.writeto("calib_mat/M_DM0_2_M_DM1.fits", M_DM0_2_M_DM1, overwrite = True)
    
    pfits.writeto("calib_mat/V_DM0_2_V_DM1.fits", V_DM0_2_V_DM1, overwrite = True)
    
    pfits.writeto("calib_mat/V_DM1_2_V_DM0.fits", V_DM1_2_V_DM0, overwrite = True)

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())


