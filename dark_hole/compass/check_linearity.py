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

    bool_DH = False
    bool_calib = True
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

    if bool_DH:
        dark_hole_phase = pfits.getdata('dark_hole.fits')
        dark_hole_phase *= 1.65/2/np.pi
        supervisor.tel.set_input_phase(dark_hole_phase)
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

    n_points = 20
    start = -1000
    end = 1000
    x = np.linspace(start, end, num=n_points)
    M_s = np.zeros(n_points)
    M_p = np.zeros(n_points)
    mode = 0
    for i in range(n_points):
    # supervisor.rtc.set_command(0,(B[:,1]+B[:,2])*np.sqrt(np.sum(pupil))/1000*1)
        supervisor.rtc.set_command(0,M2V[:,mode]/phase2nm*x[i])
        supervisor.next()
        supervisor.next()
        supervisor.next()
        phase = supervisor.target.get_tar_phase(0, pupil = True)
        slopes = supervisor.rtc.get_slopes(0)
        M_p[i] = (P2M@phase[pupil ==1]*phase2nm)[mode]
        M_s[i] = (S2M@slopes*phase2nm)[mode]

    pupil = supervisor.get_s_pupil()
    opd = np.std(phase, where = pupil ==1)
    print(opd)
    # plt.figure()
    # plt.imshow(phase)
    # plt.colorbar()
    amp = 50
    n_modes_applied = 100
    np.random.seed(2)
    command = np.zeros(n_modes)

    # command[:n_modes_applied] = (np.random.rand(n_modes_applied))*amp
    command[:n_modes_applied] = (np.random.rand(n_modes_applied)-0.5)*amp
    supervisor.rtc.set_command(0,M2V@command/phase2nm)
    supervisor.next()
    supervisor.next()
    supervisor.next()
    phase = supervisor.target.get_tar_phase(0, pupil = True)
    opd = np.std(phase, where = pupil ==1)
    print(opd)
    print(supervisor.target.get_strehl(0)[0])
    slopes = supervisor.rtc.get_slopes(0)
    phase = (P2M@phase[pupil ==1]*phase2nm)
    slopes = (S2M@slopes*phase2nm)



    plt.figure()
    plt.plot(x,M_p)
    plt.plot(x,M_s)
    plt.xlabel("tilt applied [nm]")
    plt.ylabel("tilt mesured [nm]")
    plt.legend(['measured with phase','measured with slopes'])
    # plt.title("3 lambda modulation, 20 points")
    plt.grid()
    # plt.title("no modulation")

    plt.figure()
    plt.plot(command[:n_modes_applied])
    plt.plot(phase[:n_modes_applied])
    plt.plot(slopes[:n_modes_applied])
    plt.xlabel("mode")
    plt.ylabel("amplitude [nm]")

    plt.legend(['applied','measured with phase','measured with slopes'])
    plt.grid()
    # plt.title("no modulation")


    plt.figure()
    plt.imshow(supervisor.target.get_tar_phase(0, pupil = True))
    plt.colorbar()