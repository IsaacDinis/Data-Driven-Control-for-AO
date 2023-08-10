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
        n_iter = 5000


    supervisor = Supervisor(config)
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(False) 
    phase_tilt = pfits.getdata("../data3/dist_tilt.fits")
    print(np.max(np.abs(phase_tilt)))


    pupil_diam = supervisor.config.p_geom.get_pupdiam()
    pupil_grid = make_pupil_grid(pupil_diam)
    zernike_basis = make_zernike_basis(3, 1, pupil_grid)
    tilt = zernike_basis[2].shaped
    pupil_valid = zernike_basis[0].shaped

    n_modes_DM0 = 2
    n_modes_DM1 = 2

    a = np.array([1.,-1]) 
    b = np.array([0.5,0])


    # Load command and influence matrix
    S2M_DM0 = pfits.getdata('calib_mat/S2M_DM0.fits')
    S2M_DM1 = pfits.getdata('calib_mat/S2M_DM1.fits')

    M2V_DM0 = pfits.getdata('calib_mat/M2V_DM0.fits')
    M2V_DM1 = pfits.getdata('calib_mat/M2V_DM1.fits')

    V_DM0_2_V_DM1 = pfits.getdata('calib_mat/V_DM0_2_V_DM1.fits')

    res_DM0 = np.zeros(n_iter)
    res_DM1 = np.zeros(n_iter)
    res_tilt = np.zeros(n_iter)

    command_DM0 = np.zeros(n_iter)
    command_DM1 = np.zeros(n_iter)


    #------------------------------------
    # control tilt mode
    #------------------------------------

    # res_array = np.empty((n_iter,S2M.shape[0]))
    # single_mode_res = np.empty(n_iter)

    state_mat_DM0 = np.zeros((2,2,n_modes_DM0))
    state_mat_DM1 = np.zeros((2,2,n_modes_DM1))


    bool_DMO = True
    rms_stroke = 0;


    m_pupil_size = supervisor.get_m_pupil().shape[0]
    pupil_tel_grid = make_pupil_grid(m_pupil_size)
    zernike_tel = make_zernike_basis(3, 1, pupil_tel_grid)
    tilt_tel = zernike_tel[2].shaped
    pup_valid_tel = zernike_tel[0].shaped
    tilt_record = np.dstack([tilt_tel]*5000)
    supervisor.tel.set_input_phase(tilt_record*phase_tilt[:5000])
    print(np.max(np.abs(tilt_record*phase_tilt[:5000])))

    for i in range(n_iter):
        
        slopes = supervisor.rtc.get_slopes(0)

        # modes_DM0 = np.dot(S2M_DM0,slopes)
        modes_DM1 = np.dot(S2M_DM1,slopes)
        # modes_DM1[0:n_modes_DM0] = 0

        # if  i%4==0:
        modes_DM0 = np.dot(S2M_DM0,slopes)
        state_mat_DM0[1:,:,:] = state_mat_DM0[0:-1,:,:]
        state_mat_DM0[0,0,:] = modes_DM0[0:n_modes_DM0]
        state_mat_DM0[0,1,:] = 0
        command_int_DM0 = np.dot(b,state_mat_DM0[:,0,:]) - np.dot(a,state_mat_DM0[:,1,:])
        state_mat_DM0[0,1,:] = command_int_DM0
        voltage_DM0 = -M2V_DM0[:,0:n_modes_DM0] @ command_int_DM0

        if  i%4==0:
            voltage_DM0_applied = voltage_DM0

        state_mat_DM1[1:,:,:] = state_mat_DM1[0:-1,:,:]
        state_mat_DM1[0,0,:] = modes_DM1[0:n_modes_DM1]
        state_mat_DM1[0,1,:] = 0
        command_int_DM1 = np.dot(b,state_mat_DM1[:,0,:]) - np.dot(a,state_mat_DM1[:,1,:]) 
        state_mat_DM1[0,1,:] = command_int_DM1
        voltage_DM1 = -M2V_DM1[:,0:n_modes_DM1] @ command_int_DM1
        voltage_DM1 *= 0

        if bool_DMO:
            # voltage_DM1 -= np.dot(V_DM0_2_V_DM1,voltage_DM0_applied)
            voltage = np.concatenate((voltage_DM0_applied, voltage_DM1), axis=0)
        else:
            voltage = np.concatenate((np.zeros(M2V_DM0.shape[0]), voltage_DM1), axis=0)

        supervisor.rtc.set_command(0, voltage)

        strehl = supervisor.target.get_strehl(0)

        rms_stroke += np.std(voltage_DM1)

        if i%100==0 and i > 200:
            print('s.e = {:.5f} l.e = {:.5f} \n'.format(strehl[0], strehl[1]))
        if  i%4==0:
            res_DM0[i] = modes_DM0[0]
            res_DM1[i] = modes_DM1[0]

            command_DM0[i] = voltage_DM0[0]
            command_DM1[i] = voltage_DM1[0]
        
            target_phase = supervisor.target.get_tar_phase(0,pupil=True)
            res_tilt[i] = np.sum(np.multiply(target_phase,tilt))/np.sum(pupil_valid)
        else : 
            res_DM0[i] = res_DM0[i-1]
            res_DM1[i] = res_DM1[i-1]
        

            res_tilt[i] = res_tilt[i-1]

        supervisor.next()

    rms_stroke /= n_iter
    print('rms_stroke = {:.5f} \n'.format(rms_stroke))
    start = 1000
    pfits.writeto("../data_parallel/res_DM0_alone_1kHz.fits", res_DM0[start:], overwrite = True)
    pfits.writeto("../data_parallel/res_DM0_proj_1kHz.fits", res_DM1[start:], overwrite = True)
    pfits.writeto("../data_parallel/res_tilt_DM0_1kHz.fits", res_tilt[start:], overwrite = True)

    # pfits.writeto("../data_parallel/command_DM0.fits", command_DM0[start:], overwrite = True)
    # pfits.writeto("../data_parallel/command_DM1.fits", command_DM1, overwrite = True)

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename
        embed(basename(__file__), locals())