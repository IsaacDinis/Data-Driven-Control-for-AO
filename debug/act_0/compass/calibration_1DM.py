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
from matplotlib import pyplot as plt
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
    pupil = supervisor.get_s_pupil()


    M2V_compass, _ = supervisor.basis.compute_modes_to_volts_basis("KL2V") # or "KL2V" [nvolts ,nmodes]

    xpos0 = supervisor.config.p_dms[0]._xpos # actus positions
    ypos0 = supervisor.config.p_dms[0]._ypos
    L0 = 25  # [m]
    M2V_gendron, l = KLmodes(xpos0, ypos0, L0, True) #basis on saxo stage


    n_actus = supervisor.config.p_dms[0].get_ntotact()

    n_modes = 1190
    M2V_compass = M2V_compass[:,:n_modes]
    M2V_gendron = M2V_gendron[:,:n_modes]
 
    pupil_diam = supervisor.config.p_geom.get_pupdiam()
    M2P_compass = np.zeros((int(np.sum(pupil)), n_modes))
    M2P_gendron = np.zeros((int(np.sum(pupil)), n_modes))
    
    ampli = 0.01
    slopes = supervisor.rtc.get_slopes(0)
    M2S_compass = np.zeros((slopes.shape[0], n_modes))
    M2S_gendron = np.zeros((slopes.shape[0], n_modes))

    # M2S = np.zeros((slopes.shape[0], nmodes+2))
    supervisor.atmos.enable_atmos(False)

    print('N act LODM = {:d} \n'.format(n_actus))

    print('N slopes = {:d} \n'.format(slopes.shape[0]))
    
    #-----------------------------------------------
    # compute the command matrix [nmodes , nslopes]
    #-----------------------------------------------
    for mode in range(n_modes):
        volts = M2V_compass[:,mode]*ampli
       
        supervisor.rtc.set_perturbation_voltage(0, "tmp", volts)
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/ampli
        target_phase = (supervisor.target.get_tar_phase(0,pupil=True)/ampli)[pupil ==1]

        M2S_compass[:,mode] = slopes.copy()
        M2P_compass[:,mode] = target_phase.copy()

    for mode in range(n_modes):
        volts = M2V_gendron[:,mode]*ampli
       
        supervisor.rtc.set_perturbation_voltage(0, "tmp", volts)
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/ampli
        target_phase = (supervisor.target.get_tar_phase(0,pupil=True)/ampli)[pupil ==1]

        M2S_gendron[:,mode] = slopes.copy()
        M2P_gendron[:,mode] = target_phase.copy()

    norm = np.linalg.norm(M2P_gendron,axis=0)
    M2V_gendron /= norm
    for mode in range(n_modes):
        volts = M2V_gendron[:,mode]*ampli
       
        supervisor.rtc.set_perturbation_voltage(0, "tmp", volts)
        supervisor.next()
        supervisor.next()
        slopes = supervisor.rtc.get_slopes(0)/ampli
        target_phase = (supervisor.target.get_tar_phase(0,pupil=True)/ampli)[pupil ==1]

        M2S_gendron[:,mode] = slopes.copy()
        M2P_gendron[:,mode] = target_phase.copy()

    S2M_compass = np.linalg.pinv(M2S_compass) # [nmodes , nslopes]
    noise_prop_compass = np.diagonal(S2M_compass@S2M_compass.T)

    S2M_gendron = np.linalg.pinv(M2S_gendron) # [nmodes , nslopes]
    noise_prop_gendron = np.diagonal(S2M_gendron@S2M_gendron.T)



    V2M_compass = np.linalg.pinv(M2V_compass)
    V2M_gendron = np.linalg.pinv(M2V_gendron)

    P2M_compass = np.linalg.pinv(M2P_compass)
    P2M_gendron = np.linalg.pinv(M2P_gendron)

    # pfits.writeto('calib_mat/P2M.fits', P2M, overwrite = True)
    # pfits.writeto('calib_mat/P2M_DM1.fits', P2M_DM1, overwrite = True)
    plt.figure()
    plt.semilogy(noise_prop_compass)
    plt.semilogy(noise_prop_gendron)

    p_geom = supervisor.config.p_geom
    pos_HODM = np.array([supervisor.config.p_dms[0].get_xpos(),supervisor.config.p_dms[0].get_ypos()]).T
    plt.figure()
    plt.imshow(pupil)
    plt.scatter(pos_HODM[:,0]-p_geom.get_p1(),pos_HODM[:,1]-p_geom.get_p1(), marker='.', color="blue")

    plt.figure()
    plt.plot(np.linalg.norm(M2P_compass,axis=0))
    plt.plot(np.linalg.norm(M2P_gendron,axis=0))
    command = np.zeros(n_actus)
    act = 0
    command[act]  = ampli
    supervisor.rtc.set_perturbation_voltage(0, "tmp", command)
    supervisor.next()
    supervisor.next()
    slopes = supervisor.rtc.get_slopes(0)/ampli
    target_phase = (supervisor.target.get_tar_phase(0,pupil=True)/ampli)[pupil ==1]

    # print((M2V_compass@S2M_compass@slopes)[act])
    # print((M2V_gendron@S2M_gendron@slopes)[act])
    # print((M2V_compass@P2M_compass@target_phase)[act])
    # print((M2V_gendron@P2M_gendron@target_phase)[act])

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())



        