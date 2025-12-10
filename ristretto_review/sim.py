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
#from hcipy.field import make_pupil_grid 
#from hcipy.mode_basis import make_zernike_basis 
from matplotlib import pyplot as plt
from rich.progress import track
from datetime import datetime
import os
from vibration_generator import VibrationGenerator # for vibration
#ipython -i shesha/widgets/widget_ao.py ~/Data-Driven-Control-for-AO/2DM_study/compass/compass_param.py
#V2V = np.load('../../saxo-plus/Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/V_DM0_2_V.npy')

if __name__ == "__main__":
    arguments = docopt(__doc__)
    param_file = arguments["<parameters_filename>"]

    config = ParamConfig(param_file)

    from shesha.supervisor.compassSupervisor import CompassSupervisor as Supervisor



    supervisor = Supervisor(config)
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(True) 
    #supervisor.wfs.set_noise(0, -1) #removing noise on SH WFS (-1 = no noise)

    fs = int(np.round(1/supervisor.config.p_loop.get_ittime()))
    sim_time = 5
    n_iter = int(fs*sim_time)
    #pupil_diam = supervisor.config.p_geom.get_pupdiam()
    #pupil_grid = make_pupil_grid(pupil_diam)
    #zernike_basis = make_zernike_basis(800, 1, pupil_grid)
    #tilt = zernike_basis[2].shaped
    #pupil_valid = zernike_basis[0].shaped
    now = datetime.now()
    save_path = 'results/'+ now.strftime("%Y_%m_%d_%Hh%Mm%Ss")+'_int/'
    os.mkdir(save_path)
    n_modes = 1190

    a = np.array([1.,-0.99]) 
    b = np.array([0.8,0])

    S2M = pfits.getdata('data/calib_mat/S2M.fits')
    M2V = pfits.getdata('data/calib_mat/M2V.fits')
    P2M = pfits.getdata('data/calib_mat/P2M.fits')


    res = np.zeros(n_iter)
    res_tilt = np.zeros(n_iter)



    res_phase = np.zeros((800,n_iter))


    state_mat = np.zeros((2,2,n_modes))

    boostrap_n_iter = 3000
    rms_stroke = 0;

    bool_vib = 0

    n_actus_DM0 = supervisor.config.p_dms[0].get_ntotact()
    M2V_saxo = pfits.getdata("data/KL_basis/M2V_saxo.fits")
    vib_tt_model_c =  pfits.open("data/vibration_models/modelevib0606tiptilt_2.fits") # load tiptil vibration continuous model
    vib_KL_model_c =  pfits.open("data/vibration_models/modelevib0606KL_2.fits") # load KL vibration continuous model
    vib_tt_gen = VibrationGenerator(vib_tt_model_c, fs)
    vib_KL_gen = VibrationGenerator(vib_KL_model_c, fs)
    vibration_voltages = np.zeros(supervisor.rtc.get_command(0).shape[0]-n_actus_DM0) # vibration command, size of saxo command
    pupil = supervisor.config.p_geom.get_spupil()
    phase2nm = np.sqrt(np.sum(pupil))/1e3
    P_sparta2compass = phase2nm*np.array([
    [1, 0, 0, 0, 0],
    [0, 0, -1, 0, 0],
    [0, -1, 0, 0, 0],
    [0, 0, 0, -1, 0],
    [0, 0, 0, 0, 1]])

    n_phase = 50
    phase_cube = np.zeros((*supervisor.target.get_tar_phase(0).shape, n_phase))

    for i in track(range(boostrap_n_iter), description="bootstrap"):
        slopes = supervisor.rtc.get_slopes(0)

        # modes_DM0 = np.dot(S2M_DM0,slopes)
        modes = np.dot(S2M,slopes)
        # modes[0:n_modes_DM0] = 0


        state_mat[1:,:,:] = state_mat[0:-1,:,:]
        state_mat[0,0,:] = modes[0:n_modes]
        state_mat[0,1,:] = 0
        command_int = np.dot(b,state_mat[:,0,:]) - np.dot(a,state_mat[:,1,:])
        command_int -= np.mean(command_int)  
        state_mat[0,1,:] = command_int
        voltage = -M2V[:,0:n_modes] @ command_int
        voltage[0] = 0
        # voltage *= 0

        if bool_vib:
            vibration_voltages[-2:] = vib_tt_gen.step()
            vibration_voltages[:-2] = M2V_saxo[:,2:7]@P_sparta2compass@vib_KL_gen.step()
            supervisor.rtc.set_perturbation_voltage(0,'vib',np.concatenate((np.zeros(n_actus_DM0),vibration_voltages)))
        supervisor.rtc.set_command(0, np.concatenate((voltage,np.zeros_like(vibration_voltages))))
        strehl = supervisor.target.get_strehl(0)
        if i%100==0 and i > 200:
            print('s.e = {:.5f} l.e = {:.5f} \n'.format(strehl[0], strehl[1]))
        supervisor.next()

    supervisor.target.reset_strehl(0)

    for i in track(range(n_iter), description="long exp"):
        
        slopes = supervisor.rtc.get_slopes(0)

        # modes_DM0 = np.dot(S2M_DM0,slopes)
        modes = np.dot(S2M,slopes)
        # modes[0:n_modes_DM0] = 0


        state_mat[1:,:,:] = state_mat[0:-1,:,:]
        state_mat[0,0,:] = modes[0:n_modes]
        state_mat[0,1,:] = 0
        command_int = np.dot(b,state_mat[:,0,:]) - np.dot(a,state_mat[:,1,:])
        command_int -= np.mean(command_int)  
        state_mat[0,1,:] = command_int
        voltage = -M2V[:,0:n_modes] @ command_int
        voltage[0] = 0
        # voltage *= 0


        if bool_vib:
            vibration_voltages[-2:] = vib_tt_gen.step()
            vibration_voltages[:-2] = M2V_saxo[:,2:7]@P_sparta2compass@vib_KL_gen.step()
            supervisor.rtc.set_perturbation_voltage(0,'vib',np.concatenate((np.zeros(n_actus_DM0),vibration_voltages)))
        supervisor.rtc.set_command(0, np.concatenate((voltage,np.zeros_like(vibration_voltages))))

        strehl = supervisor.target.get_strehl(0)

        rms_stroke += np.std(voltage)

        if i%(100)==0:
            print('s.e = {:.5f} l.e = {:.5f} \n'.format(strehl[0], strehl[1]))

        res[i] = modes[0]/557.2036425356519*6
        if i%(sim_time*fs/n_phase)==0:
            target_phase = supervisor.target.get_tar_phase(0,pupil=True)
            target_phase[target_phase!=0] -= np.mean(target_phase[target_phase!=0])
            phase_cube[:,:,int(i/(fs*sim_time/n_phase))] = target_phase
        #res_tilt[i] = np.sum(np.multiply(target_phase,tilt))/np.sum(pupil_valid)
        # res_DM0_all[:,i] = modes_DM0
        # res_all[:,i] = modes[:n_modes_DM0]

        # res_phase[:,i] = np.sum(np.multiply(target_phase.flatten(),zernike_basis), axis = 1)/np.sum(pupil_valid)



        supervisor.next()

    rms_stroke /= n_iter
    print('rms_stroke = {:.5f} \n'.format(rms_stroke))
    strehl = supervisor.target.get_strehl(0)
    print('s.e = {:.5f} l.e = {:.5f} \n'.format(strehl[0], strehl[1]))
    psf = supervisor.target.get_tar_image(0, expo_type='le')
    plt.figure()
    plt.imshow(np.log10(psf))
    plt.savefig(save_path+'psf.png')
    pfits.writeto(save_path+'psf.fits', psf, overwrite = True)
    pfits.writeto(save_path+'phase_cube.fits', phase_cube, overwrite = True)
    # pfits.writeto("../data_parallel/res_DM0_proj.fits", res_DM0, overwrite = True)
    # pfits.writeto("../data_parallel/res_alone.fits", res, overwrite = True)
    # pfits.writeto("../data_parallel/res_tilt.fits", res_tilt, overwrite = True)

    # pfits.writeto("../data_parallel/res_DM0_all_openloop.fits", res_DM0_all[:,10:], overwrite = True)
    # pfits.writeto("../data_parallel/res_all_openloop.fits", res_all[:,10:], overwrite = True)
    # pfits.writeto("../data_parallel/res_phase_openloop.fits", res_phase[1:89,10:], overwrite = True)
    plt.show()
    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename
        embed(basename(__file__), locals())