"""
script to control one mode (tilt)

Usage:
  closed_loop_tilt.py <parameters_filename> [options]

with 'parameters_filename' the path to the parameters file

Options:
  -h --help          Show this help message and exit
  -i, --interactive  keep the script interactive
  -d, --devices devices      Specify the devices
"""

from shesha.config import ParamConfig
from docopt import docopt
import numpy as np
from scipy.io import savemat, loadmat
import astropy.io.fits as pfits
from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 
from matplotlib import pyplot as plt
import controller
import controller_dd
import os
from datetime import datetime
import utils
from scipy.spatial import KDTree
from rich.progress import track
import eof
#ipython -i shesha/widgets/widget_ao.py ~/Data-Driven-Control-for-AO/2DM_study/compass/compass_param.py


if __name__ == "__main__":
    arguments = docopt(__doc__)

    param_file = arguments["<parameters_filename>"]

    config = ParamConfig(param_file)

    from shesha.supervisor.compassSupervisor import CompassSupervisor as Supervisor

    if arguments["--devices"]:
        config.p_loop.set_devices([
                int(device) for device in arguments["--devices"].split(",")
        ])
    supervisor = Supervisor(config)


    bool_atm = True
    bool_datadriven = False
    
    Ts = supervisor.config.p_loop.get_ittime()
    fs = 1/Ts
    exp_time = 100
    n_iter = int(np.ceil(exp_time/Ts))
    exp_time_bootstrap = 0.3
    n_bootstrap = int(np.ceil(exp_time_bootstrap/Ts))

    now = datetime.now()
    save_path = 'results/'+ now.strftime("%Y_%m_%d_%Hh%Mm%Ss")+'/'
    os.mkdir(save_path)

    
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(bool_atm) 

    pupil = supervisor.get_s_pupil()
    pupil_diam = supervisor.config.p_geom.get_pupdiam()
    pupil_grid = make_pupil_grid(pupil_diam,1)
    zernike_basis = make_zernike_basis(800, 1, pupil_grid)
    tilt = zernike_basis[2].shaped
    pupil_valid = zernike_basis[0].shaped

    DM0_phase_size = supervisor.dms.get_dm_shape(0).shape[0]
    pupil_grid_DM0 = make_pupil_grid(DM0_phase_size,1)
    zernike_basis_DM0 = make_zernike_basis(1, 0.8, pupil_grid_DM0)
    pupil_valid_DM0 = zernike_basis_DM0[0].shaped



    n_act_DM0 = supervisor.config.p_dms[0].get_ntotact()



    cross_act_DM0 = supervisor.config.p_dms[0].get_nact()+2
    # cross_act_DM0 = 48

    pos_HODM = np.array([supervisor.config.p_dms[0].get_xpos(),supervisor.config.p_dms[0].get_ypos()]).T


    n_modes_DM0 = 800

    a = np.array([1,-0.99]) 
    b = np.array([0.7,0])


    # Load command and influence matrix


    S2M_DM0 = pfits.getdata('calib_mat/S2M_DM0.fits')
    P2M_DM0 = pfits.getdata('calib_mat/P2M_DM0.fits')
    M2V_DM0 = pfits.getdata('calib_mat/M2V_DM0.fits')
    V2M = np.linalg.pinv(M2V_DM0)
    if bool_datadriven:
        IIR_filter = pfits.getdata('calib_mat/Kdd_matrix.fits')
    # P2M_DM0 = pfits.getdata('calib_mat/P2M_DM0.fits')


    res_DM0 = np.zeros(n_iter)
    res_tilt = np.zeros(n_iter)

    res_phase = np.zeros((800,n_iter))


    DM0_K = controller.K(1,a,b,S2M_DM0,M2V_DM0,stroke = np.inf)

    

    # res_array = np.empty((n_iter,S2M.shape[0]))
    # single_mode_res = np.empty(n_iter)


    state_mat_DM0 = np.zeros((2,2,n_modes_DM0))



    max_voltage = 2.20

 
    refresh_rate = 400
    DM0_plot = utils.phase_plot("tweeter phase", refresh_rate)

    target_plot = utils.phase_plot("target phase",refresh_rate)
    wfs_image_plot = utils.phase_plot("wfs image", refresh_rate)
    atm_plot = utils.phase_plot("atm phase", refresh_rate)

    zernike_saxo_plot = utils.zernike_plot("zernike res", refresh_rate, 200, pupil_diam,pupil, n_iter)
    modal_DM0_plot = utils.modal_plot("tweeter modal res", refresh_rate, n_modes_DM0, n_iter)
    modal_command_DM0_plot = utils.modal_plot("tweeter modal command", refresh_rate, n_modes_DM0, n_iter)



    DM0_stroke_plot = utils.DM_stroke_plot("tweeter stroke", refresh_rate, n_act_DM0, n_iter,pos_HODM,cross_act_DM0)
    DM0_deformation_plot = utils.deformation_plot("tweeter phase stroke", refresh_rate, n_iter)


    plt.ion()
    plt.show()


    
    cube_phase_framerate = 10
    cube_phase_count = 0
    cube_phase_target = np.zeros((pupil_diam,pupil_diam,int(np.ceil(exp_time/cube_phase_framerate/Ts))))
    DM0_phase = supervisor.dms.get_dm_shape(0)
    DM0_phase_shape = DM0_phase.shape
    cube_phase_HODM = np.zeros((DM0_phase_shape[0],DM0_phase_shape[1],int(np.ceil(exp_time/cube_phase_framerate/Ts))))
    rms_stroke = 0
    order = 2
    l = 20000
    K_eof = eof.eof(2,S2M_DM0, M2V_DM0,l)

    for i in range(n_bootstrap):
        slopes = supervisor.rtc.get_slopes(0)
        voltage = DM0_K.update_command(slopes)
        supervisor.rtc.set_command(0,voltage)
        supervisor.next()

    supervisor.target.reset_strehl(0)
    supervisor.corono.reset()
    error_rms = 0

    tilt_res = np.zeros(n_iter)
    tilt_applied = np.zeros(n_iter)
    phase_array = np.zeros(n_iter)
    slopes_array = np.zeros(n_iter)
    command_array = np.zeros(n_iter)
    for i in track(range(n_iter), description="long exposure"):

        slopes = supervisor.rtc.get_slopes(0)
        phase = supervisor.target.get_tar_phase(0)[pupil==1]
        phase -= np.mean(phase) 
        phase_array[i] = (P2M_DM0@phase)[0]
        slopes_array[i]= (S2M_DM0@slopes)[0]
        # voltage = DM0_K.update_command(slopes)
        if K_eof.is_trained == 0:
            voltage = DM0_K.update_command(slopes)
            K_eof.train(slopes,voltage)
        else: 
            # K_eof.save()
            break
            voltage_int = DM0_K.update_command(slopes)
            voltage_eof = K_eof.update_command(slopes)
            modal_command_int = V2M@voltage_int
            modal_command_eof = V2M@voltage_eof
        #     modal_command = np.copy(modal_command_int)
        #     modal_command[:2] = modal_command_eof[:2]
            # voltage = M2V_DM0@voltage_eof
            voltage = voltage_eof
        tilt_res[i] = (S2M_DM0@slopes)[0]
        tilt_applied[i] = (V2M@voltage)[0]
        command_array[i] = (V2M@voltage)[0]
        DM0_phase = supervisor.dms.get_dm_shape(0)

        supervisor.rtc.set_command(0,voltage)

        strehl = supervisor.target.get_strehl(0)

        rms_stroke += np.std(voltage)

        target_phase = supervisor.target.get_tar_phase(0,pupil=True)
        pupil = supervisor.get_i_pupil()

        phase_size = target_phase.shape[0]
        pupil_size = pupil.shape[0]
        dummy = int((pupil_size-phase_size)/2)
        pupil = pupil[dummy:-dummy,dummy:-dummy]
        target_phase *= pupil
        error_rms += np.std(target_phase*1e3,where = pupil==1)


        res_tilt[i] = np.sum(np.multiply(target_phase,tilt))/np.sum(pupil_valid)


            # print('s.e = {:.5f} l.e = {:.5f} \n'.format(strehl[0], strehl[1]))
            # print('error rms = {:.5f} \n'.format(error_rms/(i+1)))


        DM0_phase = supervisor.dms.get_dm_shape(0)#*pupil_valid_DM0
        # DM2_phase = supervisor.dms.get_dm_shape(2)#*pupil_valid_DM0
        atm_phase = supervisor.atmos.get_atmos_layer(0)
        target_phase = supervisor.target.get_tar_phase(0,pupil=True)
        target_phase[target_phase!=0] -= np.mean(target_phase[target_phase!=0])

        wfs_image = supervisor.wfs.get_wfs_image(0)



        DM0_plot.plot(DM0_phase)
        wfs_image_plot.plot(wfs_image)
        target_plot.plot(target_phase,'s.e = {:.5f} l.e = {:.5f} \n OPD rms = {:.5f} nm'.format(strehl[0], strehl[1], error_rms/(i+1)))
        
        
        atm_phase_size = atm_phase.shape[0]

        dummy = int((atm_phase_size-phase_size)/2)
        atm_phase = atm_phase[dummy:-dummy,dummy:-dummy]
        atm_phase *= -pupil
        atm_phase[atm_phase!=0] -= np.mean(atm_phase[atm_phase!=0])
        atm_plot.plot(atm_phase)
        zernike_saxo_plot.plot(target_phase,i)


        modal_DM0_plot.plot(DM0_K.res[:n_modes_DM0],i)
        modal_command_DM0_plot.plot(DM0_K.state_mat[0,1,:],i)

        DM0_deformation_plot.plot(np.max(DM0_phase)-np.min(DM0_phase),i)
        if i%cube_phase_framerate == 0:
            cube_phase_target[:,:,cube_phase_count] = target_phase
            cube_phase_HODM[:,:,cube_phase_count] = DM0_phase
            cube_phase_count += 1
        supervisor.next()


    modal_DM0_plot.compute_res_psd(fs)
    zernike_saxo_plot.compute_res_psd(fs)

    mode = 0
    modal_DM0_plot.plot_psd(mode,'HODM KL {:d} res'.format(mode),save_path+'HODM_{:d}_psd.png'.format(mode))
    zernike_saxo_plot.plot_psd(mode,'Zernike {:d} psd'.format(mode),save_path+'Zernike_{:d}_psd.png'.format(mode))

    mode = 1
    modal_DM0_plot.plot_psd(mode,'HODM KL {:d} res'.format(mode),save_path+'HODM_{:d}_psd.png'.format(mode))
    zernike_saxo_plot.plot_psd(mode,'Zernike {:d} psd'.format(mode),save_path+'Zernike_{:d}_psd.png'.format(mode))


    rms_stroke /= n_iter
    print('rms_stroke = {:.5f} \n'.format(rms_stroke))
    print('s.e = {:.5f} l.e = {:.5f} \n'.format(strehl[0], strehl[1]))
    print('error rms = {:.5f} \n'.format(error_rms/(i+1)))
    zernike_saxo_plot.save(save_path+'zernike_res.fits')


    modal_DM0_plot.save(save_path+'HODM_res.fits')
    modal_command_DM0_plot.save(save_path+'command.fits')
    utils.save_perf(save_path,exp_time,strehl[1],error_rms/(i+1))

    modal_DM0_plot.save_std_plot(save_path+'HODM_res_std.png')
    zernike_saxo_plot.save_std_plot(save_path+'zernike_std.png')

    modal_DM0_plot.save_psd(save_path+"HODM_res_psd.fits")
    zernike_saxo_plot.save_psd(save_path+"zernike_res_psd.fits")

    DM0_stroke_plot.save(save_path+'HODM_stroke.fits')
    DM0_stroke_plot.save_plot(save_path+'HODM_stroke.png')

    # hump_plot.save_plot(save_path+'hump_plot.png')
    DM0_deformation_plot.save_plot(save_path+'HODM_max_deformation.png')

    coroimg = supervisor.corono.get_image(coro_index=0) \
        / np.max(supervisor.corono.get_psf(coro_index=0))
    coroimgSampl = 1./supervisor.config.p_coronos[0].get_image_sampling()
    pfits.writeto(save_path+'corono.fits', coroimg, overwrite = True)

    contrastX, contrastAvg, contrastSigma, contrastMin, contrastMax = supervisor.corono.get_contrast(coro_index=0)
    plt.figure()
    plt.plot(contrastX, contrastAvg, '-o', label='intensity')
    plt.plot(contrastX, contrastSigma, '--', label='sigma ')
    plt.xlabel(r'r [$\lambda$ / D]')
    plt.title('Contrast curves after perfect coronograph')
    plt.ylabel('contrast')
    plt.yscale('log')
    plt.ylim(5e-7, 1e-2)
    plt.legend()
    plt.grid()
    plt.savefig(save_path+'contrast.png')

    pfits.writeto(save_path+'contrast.fits',np.array([contrastX,
           contrastAvg,
           contrastSigma,
           contrastMin,
           contrastMax]), overwrite = True)

    plt.figure()
    nimg = np.shape(coroimg)[0]
    coroX = np.arange(-nimg//2 * coroimgSampl,
                      nimg//2 * coroimgSampl, coroimgSampl)
    im=plt.pcolormesh(coroX, coroX, np.log10(coroimg))
    plt.colorbar(im)
    plt.xlabel(r'x [$\lambda$ / D]')   
    plt.ylabel(r'y [$\lambda$ / D]')
    plt.title('corono image')
    plt.savefig(save_path+'corono.png')
    pfits.writeto(save_path+'coroX.fits', coroX, overwrite = True)

    psf = supervisor.target.get_tar_image(0, expo_type='le')
    plt.figure()
    plt.imshow(np.log10(psf))
    plt.savefig(save_path+'psf.png')
    pfits.writeto(save_path+'psf.fits', psf, overwrite = True)

    pfits.writeto(save_path+'target_phase_cube.fits', cube_phase_target, overwrite = True)
    pfits.writeto(save_path+'HODM_phase_cube.fits', cube_phase_HODM, overwrite = True)

    mean_phase_res = np.mean(cube_phase_target,axis = 2)

    max_stroke_HODM = np.max(cube_phase_HODM,axis = 2)
    min_stroke_HODM = np.min(cube_phase_HODM,axis = 2)
    total_stroke_HODM = max_stroke_HODM-min_stroke_HODM

    pfits.writeto(save_path+'max_stroke_HODM.fits', total_stroke_HODM, overwrite = True)
    pfits.writeto(save_path+'upper_stroke_HODM.fits', max_stroke_HODM, overwrite = True)
    pfits.writeto(save_path+'lower_stroke_HODM.fits', min_stroke_HODM, overwrite = True)
    pfits.writeto(save_path+'mean_target_phase_res.fits', mean_phase_res, overwrite = True)

    plt.figure()
    plt.imshow(max_stroke_HODM)
    cbar = plt.colorbar()
    cbar.set_label(label="[um]", size=12)
    plt.title('upper stroke')
    plt.savefig(save_path+'upper_stroke_HODM.png')

    plt.figure()
    plt.imshow(min_stroke_HODM)
    cbar = plt.colorbar()
    cbar.set_label(label="[um]", size=12)
    plt.title('lower stroke')
    plt.savefig(save_path+'lower_stroke_HODM.png')

    plt.figure()
    plt.imshow(total_stroke_HODM)
    cbar = plt.colorbar()
    cbar.set_label(label="[um]", size=12)
    plt.title('maximum stroke')
    plt.savefig(save_path+'max_stroke_HODM.png')

    plt.figure()
    plt.imshow(mean_phase_res)
    cbar = plt.colorbar()
    cbar.set_label(label="[um]", size=12)
    plt.savefig(save_path+'mean_target_phase_res.png')

    pfits.writeto(save_path+'tilt_res.fits', tilt_res, overwrite = True)
    pfits.writeto(save_path+'tilt_applied.fits', tilt_applied, overwrite = True)
    pfits.writeto(save_path+'phase_array.fits', phase_array, overwrite = True)
    pfits.writeto(save_path+'command_array.fits', command_array, overwrite = True)
    pfits.writeto(save_path+'slopes_array.fits', slopes_array, overwrite = True)

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename
        embed(basename(__file__), locals())
