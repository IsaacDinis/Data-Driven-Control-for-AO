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
"""

from shesha.config import ParamConfig
from docopt import docopt
import numpy as np
from scipy.io import savemat, loadmat
import astropy.io.fits as pfits
from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 
from scipy.spatial import KDTree
from matplotlib import pyplot as plt
#ipython -i shesha/widgets/widget_ao.py ~/Data-Driven-Control-for-AO/2DM_study/compass/compass_param.py
#V2V = np.load('../../saxo-plus/Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/V_DM0_2_V_DM1.npy')

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
        n_iter = 40000


    supervisor = Supervisor(config)
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(True) 

    pupil_diam = supervisor.config.p_geom.get_pupdiam()
    pupil_grid = make_pupil_grid(pupil_diam)
    zernike_basis = make_zernike_basis(3, 1, pupil_grid)
    tilt = zernike_basis[1].shaped
    pupil_valid = zernike_basis[0].shaped

    n_modes_DM0 = 88
    n_modes_DM1 = 800

    n_actus_DM0 = supervisor.config.p_dms[0].get_ntotact()
    n_actus_DM1 = supervisor.config.p_dms[1].get_ntotact()

    a = np.array([1.,-1]) 
    b = np.array([0.5,0])


    # Load command and influence matrix
    S2M_DM0 = np.load('calib_mat/S2M_DM0.npy')
    S2M_DM1 = np.load('calib_mat/S2M_DM1.npy')

    M2V_DM0 = np.load('calib_mat/M2V_DM0.npy')
    M2V_DM1 = np.load('calib_mat/M2V_DM1.npy')

    V_DM0_2_V_DM1 = np.load('calib_mat/V_DM0_2_V_DM1.npy')

    res_DM0 = np.zeros(n_iter)
    res_DM1 = np.zeros(n_iter)
    res_tilt = np.zeros(n_iter)

    #------------------------------------
    # control tilt mode
    #------------------------------------

    # res_array = np.empty((n_iter,S2M.shape[0]))
    # single_mode_res = np.empty(n_iter)

    state_mat_DM0 = np.zeros((2,2,n_modes_DM0))
    state_mat_DM1 = np.zeros((2,2,n_modes_DM1))

    bool_dead_act = True
    rms_stroke = 0;


    HODM_act = 700
    pos_LODM = np.array([supervisor.config.p_dms[0].get_xpos(),supervisor.config.p_dms[0].get_ypos()]).T
    pos_HODM = np.array([supervisor.config.p_dms[1].get_xpos(),supervisor.config.p_dms[1].get_ypos()]).T
    kd_tree_LODM = KDTree(pos_LODM)


    d, i = kd_tree_LODM.query(pos_HODM[HODM_act,:], k=4)
    w = 1/d
    w /= np.sum(w)

    command_LODM = np.zeros(n_actus_DM0)
    command_dead_act = np.zeros(n_actus_DM0 + n_actus_DM1)
    for act in range(4):
        command_LODM[i[act]] = w[act]/0.49920276
        command_HODM = -V_DM0_2_V_DM1@command_LODM
        command_HODM[command_HODM>-0.0001] = 0
        command_dead_act += np.concatenate([command_LODM,command_HODM])
        command_LODM *= 0
    command_dead_act[n_actus_DM0+HODM_act] = 0

    for i in range(n_iter):

        slopes = supervisor.rtc.get_slopes(0)

        modes_DM1 = np.dot(S2M_DM1,slopes)


        state_mat_DM1[1:,:,:] = state_mat_DM1[0:-1,:,:]
        state_mat_DM1[0,0,:] = modes_DM1[0:n_modes_DM1]
        state_mat_DM1[0,1,:] = 0
        command_int_DM1 = np.dot(b,state_mat_DM1[:,0,:]) - np.dot(a,state_mat_DM1[:,1,:])
        command_int_DM1 -= np.mean(command_int_DM1)  
        state_mat_DM1[0,1,:] = command_int_DM1
        voltage_DM1 = -M2V_DM1[:,0:n_modes_DM1] @ command_int_DM1

        if bool_dead_act:
            voltage_dead_act = voltage_DM1[HODM_act] 
            voltage_DM1[HODM_act] = 0
            voltage = np.concatenate((np.zeros(M2V_DM0.shape[0]), voltage_DM1), axis=0) + command_dead_act*(voltage_dead_act-voltage_DM1[HODM_act])
        else:
            voltage = np.concatenate((np.zeros(M2V_DM0.shape[0]), voltage_DM1), axis=0)
        voltage -= np.mean(voltage)
        supervisor.rtc.set_command(0, voltage)

        strehl = supervisor.target.get_strehl(0)

        # rms_stroke += np.std(voltage_DM1)
        rms_stroke += np.abs(voltage[n_actus_DM0+HODM_act+1])
        res_DM1[i] = (voltage[n_actus_DM0+HODM_act+1])

        if i%100==0 and i > 200:
            print('s.e = {:.5f} l.e = {:.5f} \n'.format(strehl[0], strehl[1]))

        # res_DM0[i] = modes_DM0[0]/557.2036425356519
        # res_DM1[i] = modes_DM1[0]/557.2036425356519

        target_phase = supervisor.target.get_tar_phase(0,pupil=True)
        res_tilt[i] = -np.sum(np.multiply(target_phase,tilt))/np.sum(pupil_valid)
        supervisor.next()
    rms_stroke /= n_iter
    # print('rms_stroke = {:.5f} \n'.format(rms_stroke))
    print('max stroke = {:.5f} \n'.format(np.max(res_DM1)))
    # pfits.writeto("../data2/res_DM0_alone.fits", res_DM0, overwrite = True)
    pfits.writeto("../data4/act_dead.fits", res_DM1, overwrite = True)
    # pfits.writeto("../data3/res_tilt_DM0.fits", res_tilt, overwrite = True)
    # psf = supervisor.target.get_tar_image(0)
    # plt.imshow(np.log(psf))
    # plt.show()
    
    phase = supervisor.target.get_tar_phase(0,pupil=True)
    phase[phase!=0] -= np.mean(phase,where = pupil_valid.astype(bool))
    print(np.sum(np.abs(phase))/np.sum([phase!=0]))
    plt.imshow(phase)
    plt.colorbar()
    plt.show()

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename
        embed(basename(__file__), locals())