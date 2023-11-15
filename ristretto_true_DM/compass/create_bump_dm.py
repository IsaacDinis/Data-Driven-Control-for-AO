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
import utils
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


    
    p_geom = supervisor.config.p_geom
    p_tel =  supervisor.config.p_tel
    p_dm =  supervisor.config.p_dms[1]
    pixsize = supervisor.config.p_geom.get_pixsize()
    diam = p_tel.diam

    pos_HODM = np.array([supervisor.config.p_dms[1].get_xpos(),supervisor.config.p_dms[1].get_ypos()]).T
    pos_bump = np.array([supervisor.config.p_dms[2].get_xpos(),supervisor.config.p_dms[2].get_ypos()]).T
    plt.imshow(pupil)
    plt.scatter(pos_HODM[:,0]-p_geom.get_p1(),pos_HODM[:,1]-p_geom.get_p1(), marker='.', color="blue")
    plt.scatter(pos_bump[:,0]-p_geom.get_p1(),pos_bump[:,1]-p_geom.get_p1(), marker='.', color="green")
    plt.scatter(pos_HODM[376:378,0]-p_geom.get_p1(),pos_HODM[376:378,1]-p_geom.get_p1(), marker='.', color="red")
    plt.scatter(pos_HODM[416:418,0]-p_geom.get_p1(),pos_HODM[416:418,1]-p_geom.get_p1(), marker='.', color="red")
    print(p_dm._xpos.shape)
    xpos0 = 0.6 * pos_HODM[376,0] + 0.4 * pos_HODM[377,0]
    ypos0 = 0.7 * pos_HODM[376,1] + 0.3 * pos_HODM[416,1]

    plt.scatter(xpos0-p_geom.get_p1(),ypos0-p_geom.get_p1(), marker='.', color="red")
   
    xpos = p_dm._xpos-p_dm._n1
    ypos = p_dm._ypos-p_dm._n1
    xpos0 -= p_dm._n1
    ypos0 -= p_dm._n1
    i1 = p_dm._i1.copy()
    j1 = p_dm._j1.copy()
    influ = p_dm._influ.copy()

    # place bump

    influ0 = p_dm._influ[:,:,0]
    influ0 = np.expand_dims(influ0, axis=2)

    i10 = xpos0+i1[376]-xpos[376]+0.4*(i1[377]-xpos[377]-i1[376]+xpos[376])
    j10 = ypos0+j1[376]-ypos[376]+0.4*(j1[416]-ypos[416]-j1[376]+ypos[376])

    # p_dm._i1
    xcenter = p_geom.cent
    ycenter = p_geom.cent

    xpos = np.append(xpos,xpos0)
    ypos = np.append(ypos,ypos0)

    i1 = np.append(i1,i10)
    j1 = np.append(j1,j10)

    influ = np.append(influ,influ0,axis = 2)

    xpos += p_dm._n1
    ypos += p_dm._n1
    # i1 = xpos - 26
    # j1 = ypos - 26
    i1 += p_dm._n1
    j1 += p_dm._n1

    file_name = 'bump.fits'
    dm_custom = utils.write_dm_custom_fits(file_name,i1,j1,influ,xpos,ypos,xcenter,ycenter,pixsize,diam)

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
