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
        n_iter = 16000


    supervisor = Supervisor(config)
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(True) 



    pupil_diam = supervisor.config.p_geom.get_pupdiam()
    pupil_grid = make_pupil_grid(pupil_diam,1.2)
    zernike_basis = make_zernike_basis(800, 1, pupil_grid)
    tilt = zernike_basis[2].shaped
    pupil_valid = zernike_basis[0].shaped

    n_modes_DM0 = 88
    n_modes_DM1 = 800

    a = np.array([1.,-1]) 
    b = np.array([0.5,0])

    # saxo_imat = pfits.getdata("../../../compass/saxoplus-reference-cases/controller_integrated/results/saxo_imat.fits")
    # saxo_S2M = np.linalg.pinv(saxo_imat) # [n_modes , n_slopes]
    # S2M = saxo_S2M
    # cmat = pfits.getdata("../../../compass/saxoplus-reference-cases/controller_integrated/results/saxo_cmat.fits")
    # M2V = cmat@saxo_imat # [nvolts ,nmodes]
    # M2V = M2V[:-2,:]
    # Load command and influence matrix
    S2M_DM0 = np.load('calib_mat/S2M_DM0.npy')
    S2M_DM1 = np.load('calib_mat/S2M_DM1.npy')
    # S2M_DM1 = S2M
    M2V_DM0 = np.load('calib_mat/M2V_DM0.npy')
    M2V_DM1 = np.load('calib_mat/M2V_DM1.npy')
    # M2V_DM1 = M2V
    V_DM0_2_V_DM1 = np.load('calib_mat/V_DM0_2_V_DM1.npy')

    res_DM0 = np.zeros(n_iter)
    res_DM1 = np.zeros(n_iter)
    res_tilt = np.zeros(n_iter)

    res_DM0_all = np.zeros((n_modes_DM0,n_iter))
    res_DM1_all = np.zeros((n_modes_DM0,n_iter))
    res_phase = np.zeros((800,n_iter))
    voltage_DM1_all = np.zeros((n_modes_DM0,n_iter))
    #------------------------------------
    # control tilt mode
    #------------------------------------

    # res_array = np.empty((n_iter,S2M.shape[0]))
    # single_mode_res = np.empty(n_iter)

    state_mat_DM0 = np.zeros((2,2,n_modes_DM0))
    state_mat_DM1 = np.zeros((2,2,n_modes_DM1))
    V_DM1_2_V_DM0 = np.linalg.pinv(V_DM0_2_V_DM1)

    bool_DMO = True
    rms_stroke = 0;

    seeing        = [0.8]  # [arcsec]
    coherenceTime = [3.]  # [ms]
    # according to observing conditions
    r0Seeing = 500e-9/np.array(seeing)*180*3600/np.pi  # [m]

    windSpeed = np.copy(supervisor.config.p_atmos.windspeed)    # [m] / [s]
    equivSpeed = np.sum(supervisor.config.p_atmos.frac*(np.abs(windSpeed))**(5./3))**(3./5)
    tau0Seeing = 0.31 * r0Seeing / equivSpeed                # [s]
    for i in range(len(seeing)):
        print("r0 and tau0 with default speed:", r0Seeing[i],"m - ", tau0Seeing[i], "s")
    nscreens = supervisor.config.p_atmos.get_nscreens()
    supervisor.atmos.set_r0(r0Seeing[i])
    speedFactor = coherenceTime[0] * 1e-3 / tau0Seeing[i] # rescaling factors for speeds
    for l in range(nscreens): # rescale wind speeds to get the correct tau0  
        supervisor.atmos.set_wind(screen_index = l, windspeed=windSpeed[l]/speedFactor)
    voltage_DM0_applied = np.zeros(M2V_DM0.shape[0])
    for i in range(n_iter):
        
        slopes = supervisor.rtc.get_slopes(0)

        # modes_DM0 = np.dot(S2M_DM0,slopes)
        modes_DM1 = np.dot(S2M_DM1,slopes)
        # modes_DM1[0:n_modes_DM0] = 0

        modes_DM0 = np.dot(S2M_DM0,slopes)
        state_mat_DM0[1:,:,:] = state_mat_DM0[0:-1,:,:]
        state_mat_DM0[0,0,:] = modes_DM0[0:n_modes_DM0]
        state_mat_DM0[0,1,:] = 0
        command_int_DM0 = np.dot(b,state_mat_DM0[:,0,:]) - np.dot(a,state_mat_DM0[:,1,:])
        # command_int -= np.mean(command_int)  
        state_mat_DM0[0,1,:] = command_int_DM0
        voltage_DM0 = -M2V_DM0[:,0:n_modes_DM0] @ command_int_DM0


        if  i%4==3:

            voltage_DM0_applied = voltage_DM0

        state_mat_DM1[1:,:,:] = state_mat_DM1[0:-1,:,:]
        state_mat_DM1[0,0,:] = modes_DM1[0:n_modes_DM1]
        state_mat_DM1[0,1,:] = 0
        command_int_DM1 = np.dot(b,state_mat_DM1[:,0,:]) - np.dot(a,state_mat_DM1[:,1,:])
        state_mat_DM1[0,1,:] = command_int_DM1
        voltage_DM1 = -M2V_DM1[:,0:n_modes_DM1] @ command_int_DM1
        # voltage_DM1 *= 0

        if bool_DMO:
            voltage_DM0_applied = V_DM1_2_V_DM0@voltage_DM1
            voltage_DM1 -= V_DM0_2_V_DM1@voltage_DM0_applied
            voltage = np.concatenate((voltage_DM0_applied, voltage_DM1), axis=0)
        else:
            voltage = np.concatenate((np.zeros(M2V_DM0.shape[0]), voltage_DM1), axis=0)

        supervisor.rtc.set_command(0, voltage)

        strehl = supervisor.target.get_strehl(0)

        rms_stroke += np.std(voltage_DM1)

        if i%100==0 and i > 200:
            print('s.e = {:.5f} l.e = {:.5f} \n'.format(strehl[0], strehl[1]))

        # res_DM0[i] = modes_DM0[0]/557.2036425356519*6
        res_DM1[i] = modes_DM1[0]/557.2036425356519*6

        target_phase = supervisor.target.get_tar_phase(0,pupil=True)
        res_tilt[i] = np.sum(np.multiply(target_phase,tilt))/np.sum(pupil_valid)

        # res_DM0_all[:,i] = modes_DM0
        # res_DM1_all[:,i] = modes_DM1[:n_modes_DM0]

        # res_phase[:,i] = np.sum(np.multiply(target_phase.flatten(),zernike_basis), axis = 1)/np.sum(pupil_valid)
        # voltage_DM1_all[:,i] = voltage_DM1[:88]


        supervisor.next()

    rms_stroke /= n_iter
    print('rms_stroke = {:.5f} \n'.format(rms_stroke))

    # pfits.writeto("../data_parallel/res_DM0_proj.fits", res_DM0, overwrite = True)
    # pfits.writeto("../data_parallel/res_DM1_alone.fits", res_DM1, overwrite = True)
    # pfits.writeto("../data_parallel/res_tilt_DM1.fits", res_tilt, overwrite = True)

    # pfits.writeto("../gendron/res_DM0_all_DM1_alone.fits", res_DM0_all[:,150:], overwrite = True)
    # pfits.writeto("../gendron/res_DM1_all_DM1_alone.fits", res_DM1_all[:,150:], overwrite = True)
    # pfits.writeto("../gendron/res_phase_DM1_alone.fits", res_phase[1:89,150:], overwrite = True)
    # pfits.writeto("../gendron/voltage_DM1_all_DM1_alone.fits", voltage_DM1_all[:,150:], overwrite = True)
    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename
        embed(basename(__file__), locals())