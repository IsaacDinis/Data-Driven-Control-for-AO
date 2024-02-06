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

    xpos_hump = np.zeros(7)
    ypos_hump = np.zeros(7)

    i1_hump = np.zeros(7)
    j1_hump = np.zeros(7)

    xpos_hump[0] = 0.6 * pos_HODM[374,0] + 0.4 * pos_HODM[375,0]
    ypos_hump[0] = 0.7 * pos_HODM[374,1] + 0.3 * pos_HODM[414,1]

    xpos_hump[1] = 0.7 * pos_HODM[413,0] + 0.3 * pos_HODM[414,0]
    ypos_hump[1] = 0.7 * pos_HODM[413,1] + 0.3 * pos_HODM[455,1]

    xpos_hump[2] = 0.5 * pos_HODM[96,0] + 0.5 * pos_HODM[97,0]
    ypos_hump[2] = pos_HODM[96,1]

    xpos_hump[3] = 0.5 * pos_HODM[226,0] + 0.5 * pos_HODM[227,0]
    ypos_hump[3] = pos_HODM[226,1] 

    xpos_hump[4] = 0.5 * pos_HODM[261,0] + 0.5 * pos_HODM[262,0]
    ypos_hump[4] = 0.5 * pos_HODM[261,1] + 0.5 * pos_HODM[297,1]

    xpos_hump[5] = 0.6 * pos_HODM[1286,0] + 0.4 * pos_HODM[1287,0]
    ypos_hump[5] = 0.7 * pos_HODM[1286,1] + 0.3 * pos_HODM[1321,1]

    xpos_hump[6] = 0.3 * pos_HODM[1350,0] + 0.7 * pos_HODM[1351,0]
    ypos_hump[6] = 0.7 * pos_HODM[1350,1] + 0.3 * pos_HODM[1382,1]

    xpos0 = 0.5 * pos_HODM[260,0] + 0.5 * pos_HODM[261,0]
    ypos0 = 0.5 * pos_HODM[297,1] + 0.5 * pos_HODM[298,1]

    # plt.scatter(xpos0-p_geom.get_p1(),ypos0-p_geom.get_p1(), marker='.', color="red")
    # plt.scatter(xpos_hump -p_geom.get_p1(),ypos_hump -p_geom.get_p1(), marker='.', color="yellow")

    xpos = p_dm._xpos-p_dm._n1
    ypos = p_dm._ypos-p_dm._n1
    xpos0 -= p_dm._n1
    ypos0 -= p_dm._n1
    xpos_hump -= p_dm._n1
    ypos_hump -= p_dm._n1

    i1 = p_dm._i1.copy()
    j1 = p_dm._j1.copy()
    influ = p_dm._influ.copy()

    # place bump

    influ0 = p_dm._influ[:,:,0]
    influ0 = np.expand_dims(influ0, axis=2)

    influ_hump = p_dm._influ[:,:,:7]
    

    i10 = xpos0+i1[374]-xpos[374]+0.4*(i1[375]-xpos[375]-i1[374]+xpos[374])
    j10 = ypos0+j1[374]-ypos[374]+0.4*(j1[414]-ypos[414]-j1[374]+ypos[374])

    # j1_hump[0] = xpos_hump[0]+i1[374]-xpos[374]+0.4*(i1[375]-xpos[375]-i1[374]+xpos[374])
    # j1_hump[0] = ypos_hump[0]+j1[374]-ypos[374]+0.4*(j1[414]-ypos[414]-j1[374]+ypos[374])

    # j1_hump[1] = xpos_hump[1]+i1[413]-xpos[413]+0.4*(i1[414]-xpos[414]-i1[413]+xpos[413])
    # j1_hump[1] = ypos_hump[1]+j1[413]-ypos[413]+0.4*(j1[455]-ypos[455]-j1[413]+ypos[413])

    # j1_hump[2] = xpos_hump[2]+i1[96]-xpos[96]+0.4*(i1[97]-xpos[97]-i1[96]+xpos[96])
    # j1_hump[2] = ypos_hump[2]+j1[96]-ypos[96]+0.4*(j1[96]-ypos[96]-j1[96]+ypos[96])

    # j1_hump[3] = xpos_hump[3]+i1[226]-xpos[226]+0.4*(i1[227]-xpos[227]-i1[226]+xpos[226])
    # j1_hump[3] = ypos_hump[3]+j1[226]-ypos[226]+0.4*(j1[226]-ypos[226]-j1[226]+ypos[226])

    # j1_hump[4] = xpos_hump[4]+i1[261]-xpos[261]+0.4*(i1[262]-xpos[262]-i1[261]+xpos[261])
    # j1_hump[4] = ypos_hump[4]+j1[261]-ypos[261]+0.4*(j1[297]-ypos[297]-j1[261]+ypos[261])

    # j1_hump[5] = xpos_hump[5]+i1[1286]-xpos[1286]+0.4*(i1[1287]-xpos[1287]-i1[1286]+xpos[1286])
    # j1_hump[5] = ypos_hump[5]+j1[1286]-ypos[1286]+0.4*(j1[1321]-ypos[1321]-j1[1286]+ypos[1286])

    # j1_hump[6] = xpos_hump[6]+i1[1350]-xpos[1350]+0.4*(i1[1351]-xpos[1351]-i1[1350]+xpos[1350])
    # j1_hump[6] = ypos_hump[6]+j1[1350]-ypos[1350]+0.4*(j1[1382]-ypos[1382]-j1[1350]+ypos[1350])

    inf_offset = 37
    i1_hump[0] = xpos_hump[0] - inf_offset
    j1_hump[0] = ypos_hump[0] - inf_offset

    i1_hump[1] = xpos_hump[1] - inf_offset
    j1_hump[1] = ypos_hump[1] - inf_offset

    i1_hump[2] = xpos_hump[2] - inf_offset
    j1_hump[2] = ypos_hump[2] - inf_offset

    i1_hump[3] = xpos_hump[3] - inf_offset
    j1_hump[3] = ypos_hump[3] - inf_offset

    i1_hump[4] = xpos_hump[4] - inf_offset
    j1_hump[4] = ypos_hump[4] - inf_offset

    i1_hump[5] = xpos_hump[5] - inf_offset
    j1_hump[5] = ypos_hump[5] - inf_offset

    i1_hump[6] = xpos_hump[6] - inf_offset
    j1_hump[6] = ypos_hump[6] - inf_offset

    # p_dm._i1
    xcenter = p_geom.cent
    ycenter = p_geom.cent

    # xpos = np.append(xpos,xpos0)
    # ypos = np.append(ypos,ypos0)

    # i1 = np.append(i1,i10)
    # j1 = np.append(j1,j10)

    xpos = np.append(xpos,xpos_hump)
    ypos = np.append(ypos,ypos_hump)

    i1 = np.append(i1,i1_hump)
    j1 = np.append(j1,j1_hump)

    # influ = np.append(influ,influ0,axis = 2)
    influ = np.append(influ,influ_hump,axis = 2)

    xpos += p_dm._n1
    ypos += p_dm._n1
    i1 += p_dm._n1
    j1 += p_dm._n1

    xpos_hump += p_dm._n1
    ypos_hump += p_dm._n1

    i1_hump += p_dm._n1
    j1_hump += p_dm._n1

    file_name = 'bump.fits'
    dm_custom = utils.write_dm_custom_fits(file_name,i1_hump,j1_hump,influ[:,:,:7],xpos_hump,ypos_hump,xcenter,ycenter,pixsize,diam)

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
