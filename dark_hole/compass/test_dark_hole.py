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
from shesha.util.slopesCovariance import KLmodes
import astropy.io.fits as pfits
from matplotlib import pyplot as plt
from scipy.spatial import KDTree
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
    supervisor.atmos.enable_atmos(False)

    pupil = supervisor.get_s_pupil()



    p_geom = supervisor.config.p_geom
 
    dark_hole_phase = pfits.getdata('dark_hole.fits')
    dark_hole_phase *= 1.65/2/np.pi
    supervisor.tel.set_input_phase(dark_hole_phase)

    supervisor.next()
    supervisor.next()
    supervisor.next()
    target_phase = supervisor.target.get_tar_phase(0,pupil=True)



    # plt.savefig(save_path+'mean_target_phase_res.png')



    plt.figure("Target phase")
    plt.imshow(target_phase)
    plt.title("Target phase")
    cbar = plt.colorbar()
    cbar.set_label(label="[um]", size=12)

    psf = supervisor.target.get_tar_image(0, expo_type='le')

    nPixel = supervisor.get_i_pupil().shape[0]# [pixel] imaging pup. support
    pixelSize = supervisor.config.p_geom._pixsize  # [m] pupil pixel size
    psfSampl = 8./pixelSize/nPixel   # [pixel^{-1}]
    psfNCenter = psf.shape[0]//2
    psfNsampl = 100
    plt.figure()
    psfX = np.arange(-psfNsampl*psfSampl,
                     psfNsampl*psfSampl, psfSampl)
    im=plt.pcolormesh(psfX, psfX,
                      np.log10(psf[psfNCenter-psfNsampl:psfNCenter+psfNsampl,
                                         psfNCenter-psfNsampl:psfNCenter+psfNsampl]),vmin = -6, cmap='inferno')
    plt.xlabel(r'x [$\lambda$ / D]'+'\nAvg PSF before coronograph')
    plt.ylabel(r'y [$\lambda$ / D]')

    plt.colorbar(im)


    psf = supervisor.target.get_tar_image(0, expo_type='le')
    plt.figure()
    plt.imshow(np.log10(psf),vmin = -6, cmap='inferno')

    # plt.savefig(save_path+'mean_target_phase_res.png')

    # print(np.unravel_index(target_phase.argmax(), target_phase.shape))
    # print(np.unravel_index(DM1_phase.argmax(), DM1_phase.shape))
    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
