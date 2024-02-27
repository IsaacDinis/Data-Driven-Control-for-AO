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
    supervisor.atmos.enable_atmos(False) 

    
    p_geom = supervisor.config.p_geom
    p_tel =  supervisor.config.p_tel
    p_dm =  supervisor.config.p_dms[1]
    pixsize = supervisor.config.p_geom.get_pixsize()
    diam = p_tel.diam

    pos_HODM = np.array([supervisor.config.p_dms[1].get_xpos(),supervisor.config.p_dms[1].get_ypos()]).T
    pos_bump = np.array([supervisor.config.p_dms[2].get_xpos(),supervisor.config.p_dms[2].get_ypos()]).T
    plt.figure()
    plt.imshow(pupil)
    plt.scatter(pos_HODM[:,0]-p_geom.get_p1(),pos_HODM[:,1]-p_geom.get_p1(), marker='.', color="blue")
    plt.scatter(pos_bump[:,0]-p_geom.get_p1(),pos_bump[:,1]-p_geom.get_p1(), marker='.', color="green")

    # plt.scatter(pos_HODM[376:378,0]-p_geom.get_p1(),pos_HODM[376:378,1]-p_geom.get_p1(), marker='.', color="red")
    # plt.scatter(pos_HODM[416:418,0]-p_geom.get_p1(),pos_HODM[416:418,1]-p_geom.get_p1(), marker='.', color="red")
    # print(p_dm._xpos.shape)

    xpos_hump = np.zeros(8)
    ypos_hump = np.zeros(8)

    i1_hump = np.zeros(8)
    j1_hump = np.zeros(8)

    xpos_hump[0] = 0.5 * pos_HODM[305,0] + 0.5 * pos_HODM[306,0]
    ypos_hump[0] = 0.5 * pos_HODM[305,1] + 0.5 * pos_HODM[250,1]

    # plt.scatter(xpos_hump[0]-p_geom.get_p1(),ypos_hump[0]-p_geom.get_p1(), marker='.', color="red")

    xpos_hump[1] = 0.7 * pos_HODM[413,0] + 0.3 * pos_HODM[414,0]
    ypos_hump[1] = 0.7 * pos_HODM[413,1] + 0.3 * pos_HODM[455,1]

    # plt.scatter(xpos_hump[1]-p_geom.get_p1(),ypos_hump[1]-p_geom.get_p1(), marker='.', color="red")

    xpos_hump[2] = 0.5 * pos_HODM[371,0] + 0.5 * pos_HODM[372,0]
    ypos_hump[2] = pos_HODM[371,1]
    # plt.scatter(xpos_hump[2]-p_geom.get_p1(),ypos_hump[2]-p_geom.get_p1(), marker='.', color="red")

    xpos_hump[3] = 0.5 * pos_HODM[214,0] + 0.5 * pos_HODM[215,0]
    ypos_hump[3] = pos_HODM[214,1] 
    # plt.scatter(xpos_hump[3]-p_geom.get_p1(),ypos_hump[3]-p_geom.get_p1(), marker='.', color="red")

    xpos_hump[4] = 0.5 * pos_HODM[590,0] + 0.5 * pos_HODM[591,0]
    ypos_hump[4] = 0.5 * pos_HODM[590,1] + 0.5 * pos_HODM[631,1]
    # plt.scatter(xpos_hump[4]-p_geom.get_p1(),ypos_hump[4]-p_geom.get_p1(), marker='.', color="red")

    xpos_hump[5] = 0.6 * pos_HODM[545,0] + 0.4 * pos_HODM[546,0]
    ypos_hump[5] = 0.7 * pos_HODM[545,1] + 0.3 * pos_HODM[590,1]
    # plt.scatter(xpos_hump[5]-p_geom.get_p1(),ypos_hump[5]-p_geom.get_p1(), marker='.', color="red")

    xpos_hump[6] = 0.3 * pos_HODM[1558,0] + 0.7 * pos_HODM[1559,0]
    ypos_hump[6] = 0.7 * pos_HODM[1558,1] + 0.3 * pos_HODM[1582,1]
    # plt.scatter(xpos_hump[6]-p_geom.get_p1(),ypos_hump[6]-p_geom.get_p1(), marker='.', color="red")

    xpos_hump[7] = 0.3 * pos_HODM[1588,0] + 0.7 * pos_HODM[1589,0]
    ypos_hump[7] = 0.7 * pos_HODM[1588,1] + 0.3 * pos_HODM[1642,1]
    # plt.scatter(xpos_hump[7]-p_geom.get_p1(),ypos_hump[7]-p_geom.get_p1(), marker='.', color="red")

    plt.scatter(xpos_hump-p_geom.get_p1(),ypos_hump-p_geom.get_p1(), marker='.', color="red")
    # xpos0 = 0.5 * pos_bump[260,0] + 0.5 * pos_bump[261,0]
    # ypos0 = 0.5 * pos_bump[297,1] + 0.5 * pos_bump[298,1]

    # # plt.scatter(xpos0-p_geom.get_p1(),ypos0-p_geom.get_p1(), marker='.', color="red")
    # # plt.scatter(xpos_hump -p_geom.get_p1(),ypos_hump -p_geom.get_p1(), marker='.', color="yellow")

    # xpos = p_dm._xpos-p_dm._n1
    # ypos = p_dm._ypos-p_dm._n1
    # xpos0 -= p_dm._n1
    # ypos0 -= p_dm._n1
    xpos_hump -= p_dm._n1
    ypos_hump -= p_dm._n1

    # i1 = p_dm._i1.copy()
    # j1 = p_dm._j1.copy()
    influ = p_dm._influ.copy()
    influ[:,:,1] = influ[:,:,0]
    influ[:,:,2] = influ[:,:,0]
    influ[:,:,3] = influ[:,:,0]
    influ[:,:,4] = influ[:,:,0]
    influ[:,:,5] = influ[:,:,0]
    influ[:,:,6] = influ[:,:,0]
    influ[:,:,7] = influ[:,:,0]
    influ[:,:,8] = influ[:,:,0]

    # # place bump

    # influ0 = p_dm._influ[:,:,0]
    # influ0 = np.expand_dims(influ0, axis=2)

    # influ_hump = p_dm._influ[:,:,:7]
    

    # i10 = xpos0+i1[374]-xpos[374]+0.4*(i1[375]-xpos[375]-i1[374]+xpos[374])
    # j10 = ypos0+j1[374]-ypos[374]+0.4*(j1[414]-ypos[414]-j1[374]+ypos[374])


    # inf_offset = 32
    # x_offset = 20
    # y_offset = 14
    x_offset = 15.3
    y_offset = 15.3
    i1_hump[0] = xpos_hump[0] - x_offset
    j1_hump[0] = ypos_hump[0] - y_offset

    i1_hump[1] = xpos_hump[1] - x_offset
    j1_hump[1] = ypos_hump[1] - y_offset

    i1_hump[2] = xpos_hump[2] - x_offset
    j1_hump[2] = ypos_hump[2] - y_offset

    i1_hump[3] = xpos_hump[3] - x_offset
    j1_hump[3] = ypos_hump[3] - y_offset

    i1_hump[4] = xpos_hump[4] - x_offset
    j1_hump[4] = ypos_hump[4] - y_offset

    i1_hump[5] = xpos_hump[5] - x_offset
    j1_hump[5] = ypos_hump[5] - y_offset

    i1_hump[6] = xpos_hump[6] - x_offset
    j1_hump[6] = ypos_hump[6] - y_offset

    i1_hump[7] = xpos_hump[7] - x_offset
    j1_hump[7] = ypos_hump[7] - y_offset

    # # p_dm._i1
    xcenter = p_geom.cent
    ycenter = p_geom.cent


    # xpos = np.append(xpos,xpos_hump)
    # ypos = np.append(ypos,ypos_hump)

    # i1 = np.append(i1,i1_hump)
    # j1 = np.append(j1,j1_hump)

    # # influ = np.append(influ,influ0,axis = 2)
    # influ = np.append(influ,influ_hump,axis = 2)

    # xpos += p_dm._n1
    # ypos += p_dm._n1
    # i1 += p_dm._n1
    # j1 += p_dm._n1

    xpos_hump += p_dm._n1
    ypos_hump += p_dm._n1

    i1_hump += p_dm._n1
    j1_hump += p_dm._n1

    file_name = 'bump_50.fits'
    dm_custom = utils.write_dm_custom_fits(file_name,i1_hump,j1_hump,influ[:,:,:8],xpos_hump,ypos_hump,xcenter,ycenter,pixsize,diam)




    # xpos = p_dm._xpos
    # ypos = p_dm._ypos
    # i1 = p_dm._i1 + p_dm._n1
    # j1 = p_dm._j1 + p_dm._n1
    # influ = p_dm._influ / p_dm.unitpervolt

    # xcenter = p_geom.cent
    # ycenter = p_geom.cent

    # pitchm=None
    # if p_dm._pitch :
    #     pitchm = p_dm._pitch*pixsize
    # file_name = 'ristretto_50.fits'
    # dm_custom = utils.write_dm_custom_fits(file_name,i1,j1,influ,xpos,ypos,xcenter,ycenter,pixsize,diam,pitchm=pitchm)
    n_act_DM0 = supervisor.config.p_dms[0].get_ntotact()

    command = supervisor.rtc.get_command(0)
    command *= 0
    # command[-8:] = 1
    # command[n_act_DM0+305] = 1
    # command[-1] = 1
    supervisor.rtc.set_command(0,command)
    supervisor.next()
    supervisor.next()
    supervisor.next()
    supervisor.next()
    DM1_phase = supervisor.dms.get_dm_shape(1)
    phase = supervisor.target.get_tar_phase(0)
    print(np.unravel_index(phase.argmax(), phase.shape))
    # print(np.unravel_index(DM1_phase.argmax(), DM1_phase.shape))
    hump_pos_offset = 9
    hump_pos = np.array([[101,60,63,63,55,57,63,51],[24,34,30,20,44,41,100,103]])+hump_pos_offset
    plt.figure()
    plt.imshow(phase.T*pupil)
    plt.figure()
    plt.imshow(phase.T)
    plt.figure()
    plt.imshow(pupil)

    plt.scatter(xpos_hump-p_geom.get_p1(),ypos_hump-p_geom.get_p1(), marker='.', color="red")
    plt.figure()
    plt.imshow(supervisor.dms.get_dm_shape(1))
    plt.figure()
    plt.imshow(supervisor.dms.get_dm_shape(2))

    if arguments["--interactive"]:
        from shesha.util.ipython_embed import embed
        from os.path import inf_matname
        embed(inf_matname(__file__), locals())
