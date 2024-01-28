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
import os
from datetime import datetime
import utils

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
    supervisor = Supervisor(config)

    Ts = supervisor.config.p_loop.get_ittime()
    fs = 1/Ts
    exp_time = 1
    n_iter = int(np.ceil(exp_time/Ts))
    exp_time_bootstrap = 0.0
    n_bootstrap = int(np.ceil(exp_time_bootstrap/Ts))

    now = datetime.now()
    save_path = 'results/'+ now.strftime("%Y_%m_%d_%Hh%Mm%Ss")+'/'
    os.mkdir(save_path)

    
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(False) 

    pupil = supervisor.get_s_pupil()
    pupil_diam = supervisor.config.p_geom.get_pupdiam()
    pupil_grid = make_pupil_grid(pupil_diam,1)
    zernike_basis = make_zernike_basis(800, 1, pupil_grid)
    tilt = zernike_basis[2].shaped
    pupil_valid = zernike_basis[0].shaped

    DM1_phase_size = supervisor.dms.get_dm_shape(1).shape[0]
    pupil_grid_DM1 = make_pupil_grid(DM1_phase_size,1)
    zernike_basis_DM1 = make_zernike_basis(1, 0.8, pupil_grid_DM1)
    pupil_valid_DM1 = zernike_basis_DM1[0].shaped


    n_act_DM0 = supervisor.config.p_dms[0].get_ntotact()
    n_act_DM1 = supervisor.config.p_dms[1].get_ntotact()
    n_act_bump = supervisor.config.p_dms[2].get_ntotact()
    voltage_bump = np.zeros(n_act_bump)

    cross_act_DM0 = supervisor.config.p_dms[0].get_nact()+2
    # cross_act_DM1 = supervisor.config.p_dms[1].get_nact()+2
    cross_act_DM1 = 43
    pos_LODM = np.array([supervisor.config.p_dms[0].get_xpos(),supervisor.config.p_dms[0].get_ypos()]).T
    pos_HODM = np.array([supervisor.config.p_dms[1].get_xpos(),supervisor.config.p_dms[1].get_ypos()]).T

    n_modes_DM0 = 80
    n_modes_DM1 = 1000

    a = np.array([1.,-1]) 
    b = np.array([0.5,0])


    # Load command and influence matrix

    S2M_DM0 = pfits.getdata('calib_mat/S2M_DM0.fits')
    S2M_DM1 = pfits.getdata('calib_mat/S2M_DM1.fits')
    M2V_DM0 = pfits.getdata('calib_mat/M2V_DM0.fits')
    M2V_DM1 = pfits.getdata('calib_mat/M2V_DM1.fits')
    # P2M_DM1 = pfits.getdata('calib_mat/P2M_DM1.fits')
    # P2M_DM0 = pfits.getdata('calib_mat/P2M_DM0.fits')

    V_DM0_2_V_DM1 = pfits.getdata('calib_mat/V_DM0_2_V_DM1.fits')
    M2S_DM0 = np.linalg.pinv(S2M_DM0)
    V_DM1_2_V_DM0 = np.linalg.pinv(V_DM0_2_V_DM1)

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
    # DM1_K = controller.K(1,a,b,S2M_DM1,M2V_DM1,V_DM1_2_V_DM0,stroke = np.inf, offload_ratio = 4)
    DM0_K = controller.K(1,a,b,S2M_DM0,M2V_DM0)
    DM1_K = controller.K(1,a,b,S2M_DM1,M2V_DM1,stroke = np.inf)


    # res_array = np.empty((n_iter,S2M.shape[0]))
    # single_mode_res = np.empty(n_iter)

    state_mat_DM0 = np.zeros((2,2,n_modes_DM0))
    state_mat_DM1 = np.zeros((2,2,n_modes_DM1))


    bool_DMO = False
    bool_hump = False
    n_hump = 7
    hump_offset_HODM = 22
    hump_offset_hump_DM = 8
    hump_pos = np.array([[279,267,315,304,292,279,225],[141,152,46,91,109,392,415]])
    rms_stroke = 0;
    hump_amp = np.zeros(n_hump)

    voltage_DM0_applied = np.zeros(M2V_DM0.shape[0])
    refresh_rate = 100
    DM1_plot = utils.phase_plot("tweeter phase", refresh_rate)
    DM0_plot = utils.phase_plot("woofer phase", refresh_rate)
    target_plot = utils.phase_plot("target phase",refresh_rate)
    wfs_image_plot = utils.phase_plot("wfs image", refresh_rate)
    atm_plot = utils.phase_plot("atm phase", refresh_rate)

    zernike_saxo_plot = utils.zernike_plot("zernike res", refresh_rate, 200, pupil_diam,pupil, n_iter)
    modal_DM1_plot = utils.modal_plot("tweeter modal res", refresh_rate, n_modes_DM1, n_iter)
    modal_command_DM1_plot = utils.modal_plot("tweeter modal command", refresh_rate, n_modes_DM1, n_iter)
    modal_DM0_plot = utils.modal_plot("woofer modal res", refresh_rate, 80, n_iter)

    DM0_stroke_plot = utils.DM_stroke_plot("woofer stroke", refresh_rate, n_act_DM0, n_iter,pos_LODM,cross_act_DM0)
    DM1_stroke_plot = utils.DM_stroke_plot("tweeter stroke", refresh_rate, n_act_DM1, n_iter,pos_HODM,cross_act_DM1)
    DM1_deformation_plot = utils.deformation_plot("tweeter phase stroke", refresh_rate, n_iter)
    hump_deformation_plot = utils.deformation_plot("hump phase stroke", refresh_rate, n_iter)
    hump2_deformation_plot = utils.deformation_plot("hump DM2 phase stroke", refresh_rate, n_iter)
    hump_plot = utils.dummy_plot(n_hump,"hump plot",refresh_rate,n_iter)
    plt.ion()
    plt.show()

    # flat = pfits.getdata('calib_mat/flat.fits')
    # flat_record = np.dstack([flat])
    # supervisor.tel.set_input_phase(flat_record)

    for i in range(n_bootstrap):
        slopes = supervisor.rtc.get_slopes(0)
        if bool_DMO:
            voltage_DM1,voltage_DM0_applied = DM1_K.update_command(slopes)
        else:
            voltage_DM1 = DM1_K.update_command(slopes)
        voltage_DM1[386] = 2
        # voltage_DM0 = DM0_K.update_command(slopes)
        # voltage_DM0 = V_DM1_2_V_DM0@voltage_DM1
        # voltage_DM1[0] = 0
        # voltage_DM1[14] = 0
        # if  i%4==0:
        #     voltage_DM0_applied = voltage_DM0
        DM1_phase = supervisor.dms.get_dm_shape(1)
        hump_phase = supervisor.dms.get_dm_shape(2)
        # if  np.max(DM1_phase) - DM1_phase[323,162] > 1.8:
        #     voltage_bump[-1] = np.max(DM1_phase) - DM1_phase[323,162]-1.8
        if bool_hump:
            for j in range(n_hump):
                if  DM1_phase[tuple(hump_pos[:,j]+hump_offset_HODM)] < -0.1 :
                    voltage_bump[j] =  -DM1_phase[tuple(hump_pos[:,j]+hump_offset_HODM)]-0.1
                else:
                    voltage_bump[j] = 0
        else:
            voltage_bump *= 0

        for j in range(n_hump):
            hump_amp[j] = DM1_phase[tuple(hump_pos[:,j]+hump_offset_HODM)]+hump_phase[tuple(hump_pos[:,j]+hump_offset_hump_DM)]

        if bool_DMO:
            # voltage_DM1 -= V_DM0_2_V_DM1@voltage_DM0_applied
            voltage = np.concatenate((voltage_DM0_applied, voltage_DM1,voltage_bump), axis=0)
        else:
            voltage = np.concatenate((np.zeros(M2V_DM0.shape[0]), voltage_DM1,voltage_bump), axis=0)
        supervisor.rtc.set_command(0, voltage)
        supervisor.next()

    supervisor.target.reset_strehl(0)
    supervisor.target.reset_tar_phase(0)
    supervisor.corono.reset()
    error_rms = 0
    
    for i in range(n_iter):
        
        slopes = supervisor.rtc.get_slopes(0)
        # phase  = np.ndarray.flatten(supervisor.target.get_tar_phase(0,pupil=True))

        # voltage_DM0 = DM0_K.update_command(slopes)

        if bool_DMO:
            voltage_DM1,voltage_DM0_applied = DM1_K.update_command(slopes)
        else:
            voltage_DM1 = DM1_K.update_command(slopes)
        voltage_DM1 *= 0
        voltage_DM1[386] = 2 
        # voltage_DM1[0] = 0
        # voltage_DM1[14] = 0
        # voltage_DM0 = V_DM1_2_V_DM0@voltage_DM1

        # if  i%4==0:
        #     voltage_DM0_applied = voltage_DM0

        DM1_phase = supervisor.dms.get_dm_shape(1)
        hump_phase = supervisor.dms.get_dm_shape(2)

        # if  DM1_phase[323,162] < -0.1 :
        #     voltage_bump[-1] =  -DM1_phase[323,162]-0.1

        # if  np.max(DM1_phase) - DM1_phase[323,162] > 1.8:
        #     voltage_bump[-1] = np.max(DM1_phase) - DM1_phase[323,162]-1.8
        # else:
        #     voltage_bump[-1] = 0
        # voltage_bump[-1] = 0
        if bool_hump:
            for j in range(n_hump):
                if  DM1_phase[tuple(hump_pos[:,j]+hump_offset_HODM)] < -0.1 :
                    voltage_bump[j] =  -DM1_phase[tuple(hump_pos[:,j]+hump_offset_HODM)]-0.1
                else:
                    voltage_bump[j] = 0
                
        else:
            voltage_bump *= 0

        # for j in range(n_hump):
        #     hump_amp[j] = DM1_phase[tuple(hump_pos[:,j]+hump_offset_HODM)]+hump_phase[tuple(hump_pos[:,j]+hump_offset_hump_DM)]


        if bool_DMO:
            # voltage_DM1 -= V_DM0_2_V_DM1@voltage_DM0_applied
            voltage = np.concatenate((voltage_DM0_applied, voltage_DM1,voltage_bump), axis=0)
        else:
            voltage = np.concatenate((np.zeros(M2V_DM0.shape[0]), voltage_DM1,voltage_bump), axis=0)

        supervisor.rtc.set_command(0, voltage)

        strehl = supervisor.target.get_strehl(0)

        rms_stroke += np.std(voltage_DM1)

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

        DM0_phase = supervisor.dms.get_dm_shape(0)
        DM1_phase = supervisor.dms.get_dm_shape(1)#*pupil_valid_DM1
        DM2_phase = supervisor.dms.get_dm_shape(2)#*pupil_valid_DM1
        atm_phase = supervisor.atmos.get_atmos_layer(0)
        target_phase = supervisor.target.get_tar_phase(0,pupil=True)
        target_phase[target_phase!=0] -= np.mean(target_phase[target_phase!=0])
        wfs_image = supervisor.wfs.get_wfs_image(0)

        # DM0_plot.plot(DM0_phase,'hump = {:.5f} '.format(target_phase[301,140]))
        # DM1_plot.plot(DM1_phase,'hump = {:.5f} '.format(np.max(DM1_phase) - DM1_phase[323,162]))
        wfs_image_plot.plot(wfs_image)
        target_plot.plot(target_phase,'s.e = {:.5f} l.e = {:.5f} \n OPD rms = {:.5f} nm'.format(strehl[0], strehl[1], error_rms/(i+1)))
        
        
        atm_phase_size = atm_phase.shape[0]

        dummy = int((atm_phase_size-phase_size)/2)
        atm_phase = atm_phase[dummy:-dummy,dummy:-dummy]
        atm_phase *= -pupil
        atm_phase[atm_phase!=0] -= np.mean(atm_phase[atm_phase!=0])
        atm_plot.plot(atm_phase)
        zernike_saxo_plot.plot(target_phase,i)

        modal_DM0_plot.plot(DM0_K.res,i)
        modal_DM1_plot.plot(DM1_K.res[:n_modes_DM1],i)
        modal_command_DM1_plot.plot(DM1_K.state_mat[0,1,:],i)

        # DM0_stroke_plot.plot(voltage_DM0_applied,i)
        # DM1_stroke_plot.plot(voltage_DM1,i)
        DM1_deformation_plot.plot(np.max(DM1_phase)-np.min(DM1_phase),i)
        # hump_deformation_plot.plot(np.max(DM1_phase)- DM1_phase[323,162],i)
        # hump_deformation_plot.plot(DM1_phase[225+22,415+22],i)
        # hump2_deformation_plot.plot(DM2_phase[225+8,415+8],i)
        # hump_plot.plot(hump_amp,i)
        supervisor.next()

    modal_DM0_plot.compute_res_psd(fs)
    modal_DM1_plot.compute_res_psd(fs)
    zernike_saxo_plot.compute_res_psd(fs)

    mode = 0
    modal_DM0_plot.plot_psd(mode,'LODM KL {:d} res'.format(mode),save_path+'LODM_{:d}_psd.png'.format(mode))
    modal_DM1_plot.plot_psd(mode,'HODM KL {:d} res'.format(mode),save_path+'HODM_{:d}_psd.png'.format(mode))
    zernike_saxo_plot.plot_psd(mode,'Zernike {:d} psd'.format(mode),save_path+'Zernike_{:d}_psd.png'.format(mode))

    mode = 1
    modal_DM0_plot.plot_psd(mode,'LODM KL {:d} res'.format(mode),save_path+'LODM_{:d}_psd.png'.format(mode))
    modal_DM1_plot.plot_psd(mode,'HODM KL {:d} res'.format(mode),save_path+'HODM_{:d}_psd.png'.format(mode))
    zernike_saxo_plot.plot_psd(mode,'Zernike {:d} psd'.format(mode),save_path+'Zernike_{:d}_psd.png'.format(mode))


    rms_stroke /= n_iter
    print('rms_stroke = {:.5f} \n'.format(rms_stroke))
    print('s.e = {:.5f} l.e = {:.5f} \n'.format(strehl[0], strehl[1]))
    print('error rms = {:.5f} \n'.format(error_rms/(i+1)))
    zernike_saxo_plot.save(save_path+'zernike_res.fits')

    modal_DM0_plot.save(save_path+'LODM_res.fits')
    modal_DM1_plot.save(save_path+'HODM_res.fits')
    utils.save_perf(save_path,exp_time,strehl[1],error_rms/(i+1))

    modal_DM1_plot.save_std_plot(save_path+'HODM_res_std.png')
    modal_DM0_plot.save_std_plot(save_path+'LODM_res_std.png')
    zernike_saxo_plot.save_std_plot(save_path+'zernike_std.png')

    modal_DM0_plot.save_psd(save_path+"LODM_res_psd.fits")
    modal_DM1_plot.save_psd(save_path+"HODM_res_psd.fits")
    zernike_saxo_plot.save_psd(save_path+"zernike_res_psd.fits")

    DM0_stroke_plot.save(save_path+'LODM_stroke.fits')
    DM1_stroke_plot.save(save_path+'HODM_stroke.fits')

    DM0_stroke_plot.save_plot(save_path+'LODM_stroke.png')
    DM1_stroke_plot.save_plot(save_path+'HODM_stroke.png')

    hump_plot.save_plot(save_path+'hump_plot.png')
    DM1_deformation_plot.save_plot(save_path+'HODM_max_deformation.png')

    coroimg = supervisor.corono.get_image(coro_index=0) \
        / np.max(supervisor.corono.get_psf(coro_index=0))
    coroimgSampl = 1./supervisor.config.p_coronos[0].get_image_sampling()
    pfits.writeto(save_path+'corono.fits', coroimg, overwrite = True)

    plt.figure()
    nimg = np.shape(coroimg)[0]
    coroX = np.arange(-nimg//2 * coroimgSampl,
                      nimg//2 * coroimgSampl, coroimgSampl)
    im=plt.pcolormesh(coroX, coroX, coroimg)
    plt.colorbar(im)
    plt.xlabel(r'x [$\lambda$ / D]')   
    plt.ylabel(r'y [$\lambda$ / D]')
    plt.title('corono image')
    plt.savefig(save_path+'corono.png')

    psf = supervisor.target.get_tar_image(0, expo_type='le')
    plt.figure()
    plt.imshow(np.log10(psf))
    plt.savefig(save_path+'psf.png')

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import basename
        embed(basename(__file__), locals())