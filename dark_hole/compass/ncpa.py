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
from matplotlib import pyplot as plt
import os
from datetime import datetime
from shesha.util.slopesCovariance import KLmodes

if __name__ == "__main__":

  
    bool_calib = False
    bool_NCPA = False
    arguments = docopt(__doc__)

    param_file = arguments["<parameters_filename>"]

    config = ParamConfig(param_file)

    from shesha.supervisor.compassSupervisor import CompassSupervisor as Supervisor

    if arguments["--devices"]:
        config.p_loop.set_devices([
                int(device) for device in arguments["--devices"].split(",")
        ])
    supervisor = Supervisor(config)
    
    supervisor.rtc.open_loop(0) # disable implemented controller
    supervisor.atmos.enable_atmos(False) 
    pupil = supervisor.config.p_geom.get_spupil()
    n_pixels = int(np.sum(pupil))
    phase2nm = 1e3/np.sqrt(np.sum(pupil))
    p_geom = supervisor.config.p_geom
    p_dm =  supervisor.config.p_dms[0]
    
    n_modes = 400
    M2V = pfits.getdata("saxoplus_M2V.fits")

    # xpos = supervisor.config.p_dms[0]._xpos # actus positions
    # ypos = supervisor.config.p_dms[0]._ypos
    # L0 = 25  # [m]
    # M2V, l = KLmodes(xpos, ypos, L0, True) #basis on saxo stage
    # M2V = M2V[:,:n_modes]

    ampli = 1.0e-2
    n_slopes = supervisor.rtc.get_slopes(0).shape[0]
    imat = np.zeros((n_slopes, n_modes))
    imatp = np.zeros((n_pixels, n_modes))


    if bool_NCPA:
        ncpa_phase = pfits.getdata('dark_hole.fits')
        ncpa_phase_s = pfits.getdata('dark_hole_s_size.fits')
        ncpa_phase *= 1.65/2/np.pi
        supervisor.wfs.set_wfs_phase(0,ncpa_phase)

        supervisor.wfs.compute_wfs_image(0, noise=False)

        supervisor.rtc.do_centroids(0)

        rs = supervisor.rtc.get_slopes(0)
        supervisor.tel.set_input_phase(ncpa_phase)
        supervisor.next()
        supervisor.next()
        supervisor.next()


    if bool_calib:
        for mode in range(n_modes):
            volts = M2V[:, mode] * ampli; 
            supervisor.rtc.set_perturbation_voltage(0, "tmp", volts)
            supervisor.next()
            supervisor.next()
            supervisor.next()
            s = supervisor.rtc.get_slopes(0)/ampli
            p = supervisor.target.get_tar_phase(0)[pupil==1]/ampli
            p -= np.mean(p) 
            
            imat[:, mode] = s.copy()
            imatp[:, mode] = p.copy()

        # P2M = np.linalg.pinv(imatp)
        # S2M = np.linalg.pinv(imat)
        P2M = np.dot(np.linalg.inv(np.dot(imatp.T, imatp)), imatp.T)
        S2M = np.dot(np.linalg.inv(np.dot(imat.T, imat)), imat.T)
        
        pfits.writeto('P2M.fits', P2M, overwrite = True)
        pfits.writeto('S2M.fits', S2M, overwrite = True)
    else:
        P2M = pfits.getdata("P2M.fits")
        S2M = pfits.getdata("S2M.fits")


    supervisor.rtc.reset_ref_slopes(0)
    supervisor.rtc.set_perturbation_voltage(0, "tmp", M2V[:, 0]*0)
    supervisor.next()
    supervisor.next()
    supervisor.next()
    s = supervisor.rtc.get_slopes(0)

    plt.figure()
    plt.imshow(supervisor.target.get_tar_phase(0, pupil = True))
    cbar = plt.colorbar()
    cbar.set_label(label="[um]", size=12)
    plt.title('phase')

    psf = supervisor.target.get_tar_image(0)
    nPixel = supervisor.get_i_pupil().shape[0]# [pixel] imaging pup. support 
    pixelSize = supervisor.config.p_geom._pixsize  # [m] pupil pixel size
    diameter = supervisor.config.p_tel.diam
    psfSampl = diameter/pixelSize/nPixel   # [pixel^{-1}]
    psfNCenter = psf.shape[0]//2
    psfNsampl = 30

    plt.figure()
    psfX = np.arange(-psfNsampl*psfSampl,
                     psfNsampl*psfSampl, psfSampl)
    im=plt.pcolormesh(psfX, psfX,
                      np.log10(psf[psfNCenter-psfNsampl:psfNCenter+psfNsampl,
                                         psfNCenter-psfNsampl:psfNCenter+psfNsampl]), vmin = -8, vmax = 0)
    plt.xlabel(r'x [$\lambda$ / D]'+'\nAvg PSF before coronograph')
    plt.ylabel(r'y [$\lambda$ / D]')

    plt.title('PSF')
    cbar = plt.colorbar()
    cbar.set_label(label="contrast", size=12)


    plt.figure()
    plt.imshow(supervisor.wfs.get_wfs_phase(0))

    supervisor.target.set_ncpa_tar(0,ncpa_phase_s)
    supervisor.next()
    supervisor.next()
    supervisor.next()

    supervisor.wfs.set_ncpa_wfs(0,ncpa_phase*0)
    supervisor.next()
    supervisor.next()
    supervisor.next()

    s = supervisor.rtc.get_slopes(0)
    supervisor.rtc.set_ref_slopes(s)
    supervisor.next()
    supervisor.next()
    supervisor.next()
    s = supervisor.rtc.get_slopes(0)