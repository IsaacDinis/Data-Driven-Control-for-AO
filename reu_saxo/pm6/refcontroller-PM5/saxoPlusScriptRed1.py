"""
Author: F.V - Modified by C. Bechet (2023.01.03)                                          

This script is intended to be used with the saxo+ Manager ONLY!

The saxo+ manager is instantiated with the reference configuration files 
saxo.py and saxo+.py

The file saxo.py contains the saxo parameter configuration file.
The file saxo+.py contains the saxo+ parameter configuration file.

COMPASS methods can afetrwards and during the script, through the saxoplus hooks be accessed using either variables:

<saxoplusmanager>, the main SAXO+ manager (which handles the synchronisation between the 2 AO stages)
<saxo>           , the First stage compass supervisor
<saxoplus>       , the Second stage compass supervisor


Usage:
  saxoPlusScript.py <controllerpath> <dataDir>

Arguments:
  <controllerpath> used to specify the path to the controller scripts for the hooks.
  <dataDir>        used to store the results of the simulations

Example:
  ipython -i saxoPlusScript.py <controllerpath> <dataDir>
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pfits
import os
import datetime
from shesha.util.slopesCovariance import KLmodes
from saxoPlusUtils import *
from saxoPlusOutputs import Outputs, OutputsType

from tqdm import trange

from shesha.config import ParamConfig
from twoStagesManager import TwoStagesManager
from stage1Supervisor import Stage1Supervisor
from stage2Supervisor import Stage2Supervisor

from docopt import docopt
arguments = docopt(__doc__)
controllerpath = str(arguments['<controllerpath>'])
dataDir = str(arguments['<dataDir>'])

"""
Simulation parameters
"""
# observing conditions
seeing        = [0.7]  # [arcsec]
coherenceTime = [3]  # [ms]
photFluxStage1 = [29200] # [photons / m2 / s]
photFluxStage2 = [858000]  # [photons / m2 / s]

# define bootstrap time
bootstrap_time = 0.1  # [s]                                                              
# define long exposure time
exposure_time = 5.  # [s]

"""
Loading of the configurations of first and second stages.
Setup of the saxoplusmanager.
USEFUL VARIABLES are: saxo, saxoplus and saxoplusmanager.
"""

config1=ParamConfig("../cases/saxoRed1.py")
config2=ParamConfig("../cases/saxo+Red1.py")

# Before going further, an integer frequency_ratio must exist between the two stages.
frequency_ratio = sternBrocotRatio(config1.p_loop.ittime, config2.p_loop.ittime, rel_tol=1e-2)

# need to rescale parameters in frame units according to new frequency: ittime, niter and delay.
newittime = config1.p_loop.ittime/frequency_ratio[0]
newniter = int(config1.p_loop.niter*frequency_ratio[0])
config1.p_loop.set_ittime(newittime)
config1.p_loop.set_niter(newniter)
config2.p_loop.set_ittime(newittime)
config2.p_loop.set_niter(newniter)

# Adapt the RON noise in both stages to account for sub-exposure times
config1.p_wfss[0].set_noise(np.sqrt(config1.p_wfss[0].get_noise()**2 / frequency_ratio[0]))                                                                  
config2.p_wfss[0].set_noise(np.sqrt(config2.p_wfss[0].get_noise()**2 / frequency_ratio[1]))


nControl = range(len(config1.p_controllers))
for ncontrol in nControl:
    config1.p_controllers[ncontrol].set_delay(config1.p_controllers[ncontrol].get_delay()*frequency_ratio[0])

nControl = range(len(config2.p_controllers))
for ncontrol in nControl:
    config2.p_controllers[ncontrol].set_delay(config2.p_controllers[ncontrol].get_delay()*frequency_ratio[1])
    
effFreqs = 1./(newittime*frequency_ratio)
print("Effectively simulated frequencies:", effFreqs)
        
first_stage = Stage1Supervisor(config1)
second_stage = Stage2Supervisor(config2)

saxoplusmanager = TwoStagesManager(first_stage, second_stage,
                                   frequency_ratio=frequency_ratio[0],
                                   frequency_ratio_denom=frequency_ratio[1])

saxo = saxoplusmanager.first_stage
saxoplus = saxoplusmanager.second_stage
saxo.manager=saxoplusmanager
saxoplus.manager=saxoplusmanager
saxoplusmanager.ittime = newittime

""" 
Add some extra required configuration on supervisors: clipping on BOSTON DM.
"""
vClip=3.2 # +/- 3.2 microns on each actuator
saxoplus.rtc._rtc.d_control[0].set_comRange(-vClip, vClip)

"""
Add the hooks
"""
hooks_available = ['Config', 'CalibOnSource','CalibOnSky',
                   'Stage1NewSlopes', 'Stage2NewSlopes']
nhooks = len(hooks_available)
hooks=[] # initialize hooks and list of hooks
hooks_list = []
for i in range(nhooks):
    if (os.path.exists(controllerpath+'saxoPlus'+hooks_available[i]+'.py')):
        hooks_list.append(hooks_available[i])
        hooks.append(controllerpath+'saxoPlus'+hooks_available[i]+'.py')
saxoplusmanager.hooks_list = hooks_list
saxoplusmanager.hooks = hooks
saxoplusmanager.controllerpath = controllerpath



"""                                                                                      
Call for any further configuration required by the controller
using saxoPlusConfig.py if necessary                                           
"""
if 'Config' in saxoplusmanager.hooks_list:
    hook_index = saxoplusmanager.hooks_list.index('Config')
    exec(open(saxoplusmanager.hooks[hook_index]).read())


"""
Prepare for On-Source Calibration.
And call saxoPlusCalibOnSource.py if necessary
"""
saxo.atmos.enable_atmos(False) # disabling turbulence
if 'CalibOnSource' in saxoplusmanager.hooks_list:
    hook_index = saxoplusmanager.hooks_list.index('CalibOnSource')
    exec(open(saxoplusmanager.hooks[hook_index]).read())


"""
Prepare for closed-loop runs
"""
saxo.atmos.enable_atmos(True) # enabling turbulence
saxo.rtc.close_loop(0) # closing loop on first stage
saxoplus.rtc.close_loop(0) # closing loop on second stage

# configure bootstrap and long exposure                                                  
bootstrap_iter = int(bootstrap_time / newittime)
exposure_iter = int(exposure_time / newittime)

# computing Fried parameter and compass 'gs_mag'
# according to observing conditions
r0Seeing = 500e-9/np.array(seeing)*180*3600/np.pi  # [m]

windSpeed = np.copy(saxo.config.p_atmos.windspeed)    # [m] / [s]
equivSpeed = np.sum(saxo.config.p_atmos.frac*(np.abs(windSpeed))**(5./3))**(3./5)
tau0Seeing = 0.31 * r0Seeing / equivSpeed                # [s]
for i in range(len(seeing)):
    print("r0 and tau0 with default speed:", r0Seeing[i],"m - ", tau0Seeing[i], "s")
nscreens = saxo.config.p_atmos.get_nscreens()

# creating strehl.txt file containing strehl
SimuDataBaseFile = dataDir+'simus.txt'
if os.path.exists(SimuDataBaseFile):
    dbfile = open(SimuDataBaseFile, 'a')
else:
    dbfile = open(SimuDataBaseFile, 'a')
    line = 'seeing [arcsec] | req. r0 [m] | eff. r0 [m] | req. tau0 [ms] | eff. tau0 [ms] |req. flux stage1 [photons/m2/s] | eff. flux [photons/subap/frame] | req. flux stage2  [photons/m2/s] | eff. flux [photons/pix/frame] | avg strehl LE SAXO | sigma strehl SE SAXO | avg strehl LE SAXO+ | sigma strehl SE SAXO+ |'
    dbfile.write(line + '\n')


# loops
for i in range(len(seeing)):
    for j in range(len(coherenceTime)):
        speedFactor = coherenceTime[j] * 1e-3 / tau0Seeing[i] # rescaling factors for speeds
        print("Speed rescaling factor:", speedFactor)

        for k in range(len(photFluxStage1)):
            for kk in range(len(photFluxStage2)):

                prefix=createDirForOutputs(dataDir)
                # setting flux observing conditions
                saxo.set_gsmag_from_photFlux(photFluxStage1[k])
                saxoplus.set_gsmag_from_photFlux(photFluxStage2[kk])

                # setting atmospheric observing conditions                      
                saxo.atmos.set_r0(r0Seeing[i])
                for l in range(nscreens): # rescale wind speeds to get the correct tau0  
                    saxo.atmos.set_wind(screen_index = l, windspeed=windSpeed[l]/speedFactor)
                    
                # RESET FOR A NEW STARTING CLOSED-LOOP SIMULATION
                saxo.atmos.reset_turbu()
                saxo.wfs.reset_image()
                saxo.dms.reset_dm()
                saxo.rtc.reset_command()
                saxoplus.wfs.reset_image()
                saxoplus.dms.reset_dm()
                saxoplus.rtc.reset_command()
                saxoplusmanager.iterations = 0 # reset iterations to zero to reset WFS integration    
                # starting bootstrap
                print("starting bootstrap")
                for e in trange(bootstrap_iter): saxoplusmanager.next()

                # starting long exposure... 
                # reset metrics arrays
                outputsKeys = ['strehl','psf','rtc','phase','corono','contrast']
                saxo.outputs = Outputs(saxo, exposure_iter, 1, outputsKeys,
                                       prefix)             
                saxoplus.outputs = Outputs(saxoplus, exposure_iter, 2,
                                           outputsKeys, prefix)

                print("starting long exposure")
                plt.ion()
                saxo.target.reset_strehl(tar_index = 0) # resetting long exp. SR on saxo
                saxoplus.target.reset_strehl(tar_index = 0) # resetting long exp. SR on saxo+
                saxo.corono.reset()  # reset coronagraph integration
                saxoplus.corono.reset()

                for e in trange(exposure_iter):
                    saxoplusmanager.next() # exposure...
                    saxo.outputs.update(saxo, e)
                    saxoplus.outputs.update(saxoplus, e)

                # at the end of the exposure...    
                saxo.outputs.finalize(saxo)
                saxoplus.outputs.finalize(saxoplus)

                # gather effectively simulated parameters
                r0eff = saxo.config.p_atmos.get_r0() # [m]                      
                tau0eff = 0.31 * r0eff / (np.sum(saxo.config.p_atmos.frac*(np.abs(saxo.config.p_atmos.get_windspeed()))**(5./3))**(3./5))            # [s]
                nphotons_sh = saxo.compute_nphotons(photFluxStage1[k],
                                                    effFreqs[0])
                nphotons_pyr = saxoplus.compute_nphotons(photFluxStage2[kk],
                                                         effFreqs[1])
                line = '{} {:.2f} {:.2f} {} {:.2f} {:.0f} {:.0f} {:.0f} {:.0f} {:.5f} {:.5f} {:.5f} {:.5f} {:s}'.format(seeing[i], r0Seeing[i], r0eff, coherenceTime[j], tau0eff*1e3, photFluxStage1[k], nphotons_sh, photFluxStage2[kk], nphotons_pyr, np.mean(saxo.outputs.strehlVal[1,:]), np.sqrt(np.var(saxo.outputs.strehlVal[0,:])), np.mean(saxoplus.outputs.strehlVal[1,:]), np.sqrt(np.var(saxoplus.outputs.strehlVal[0,:])), prefix)
                dbfile.write(line + '\n')

                saveOutputs(saxoplusmanager, prefix)
                outputstitle='r0={:.2f}m tau0={:.1f}ms freqs {:.0f} & {:.0f} Hz\n SH {:.0f}phot/sub/frame \n Pyr {:.0f}phot/pix/frame'.format(r0eff, tau0eff*1e3, effFreqs[0], effFreqs[1], nphotons_sh, nphotons_pyr)
                saxo.outputs.plot(outputstitle=outputstitle)
                saxoplus.outputs.plot(outputstitle=outputstitle)
                
                print("seeing = ", seeing[i])
                print("coherence time = ", coherenceTime[j])
                print("flux stage1 = ", photFluxStage1[k])
                print("flux stage2 = ", photFluxStage2[kk])
                print("SAXO final LE SR = ",
                      saxo.outputs.strehlVal[1,exposure_iter-1]) 
                print("SAXO+ final LE SR = " ,
                      saxoplus.outputs.strehlVal[1,exposure_iter-1])
                
dbfile.close()

