#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Isaac Dinis
"""
#ipython -i shesha/widgets/widget_ao.py ~/Data-Driven-Control-for-AO/2WFS_study/compass/param_1WFS.py

import shesha.config as conf
import numpy as np

#simul_name = "saxo+"
Ts = 1./1000.

# loop
p_loop = conf.Param_loop()
p_loop.set_niter(5000)                      # number of loops will be overwritten by SAXO
p_loop.set_ittime(Ts)     # /!\ second stage frequency to be set above


# geom
p_geom = conf.Param_geom()
p_geom.set_zenithangle(0.)      # /!\ keep it 0 until we know what it does in COMPASS

# tel
p_tel = conf.Param_tel()
p_tel.set_diam(8.0)            # /!\  VLT diameter
p_tel.set_cobs(0)           # /!\  central obstruction
p_tel.set_type_ap('VLT-NoObs')      # /!\  VLT pupil
p_tel.set_t_spiders(0.0)   # /!\  spider width = 5 cm

# atmos
p_atmos = conf.Param_atmos()  # /!\ Dummy atmosphere for setup of WFS ...
p_atmos.set_r0(0.14)          # /!\ Fried parameter will be overwritten by the SAXO.PY value
p_atmos.set_nscreens(1)       # /!\ Number of layers
p_atmos.set_frac([1.0])       # /!\ Fraction of atmosphere (100% = 1)
p_atmos.set_alt([0.0])        # /!\ Altitude(s) in meters
p_atmos.set_windspeed([10])   # /!\ Wind speed of layer(s) in m/s
p_atmos.set_winddir([45])     # /!\ Wind direction in degrees
p_atmos.set_L0([25])          # /!\ Outer scale in meters

# target
p_target = conf.Param_target()
p_targets = [p_target]
p_target.set_xpos(0.)         # /!\ On axis
p_target.set_ypos(0.)         # /!\ On axis
p_target.set_Lambda(1.65)     # /!\ H Band
p_target.set_mag(6.)          # /!\

# wfs
p_wfs0 = conf.Param_wfs(roket=True)
p_wfss = [p_wfs0]

p_wfs0.set_type("pyrhr")        # /!\ pyramid
p_wfs0.set_nxsub(16)            #     number of pixels
p_wfs0.set_fracsub(0.5)      #     threshold on illumination fraction for valid pixel
p_wfs0.set_Lambda(1.2)          #     wavelength
p_wfs0.set_gsmag(6.)
p_wfs0.set_zerop(1.e11)
p_wfs0.set_optthroughput(1)
p_wfs0.set_noise(-1)           #     readout noise
p_wfs0.set_xpos(0.)             # /!\ On axis
p_wfs0.set_ypos(0.)             # /!\ On axis
rMod = 3.                       # Modulation radius, in lam/D units
p_wfs0.set_pyr_ampl(rMod)
nbPtMod = int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.) * 4)
p_wfs0.set_pyr_npts(nbPtMod) 
p_wfs0.set_pyr_pup_sep(p_wfs0.nxsub) # separation between the 4 images of the pyramid 
p_wfs0.set_fstop("round")
p_wfs0.set_fssize(1.5)          # Size of the field stop
p_wfs0.set_atmos_seen(1)        # /!\


p_dm0 = conf.Param_dm()       # /!\
p_dm0.set_type("tt")         # /!\
p_dm0.set_alt(0.)            # /!\
p_dm0.set_unitpervolt(1.)    # /!\
p_dm0.set_push4imat(0.18)    #     to displace about half a pixel
p_dms = [p_dm0]

# centroiders
p_centroider0 = conf.Param_centroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)           # /!\
p_centroider0.set_type("maskedpix")

# controllers
p_controller0 = conf.Param_controller()
p_controllers = [p_controller0]

#p_controller0.set_type("generic")   # /?\ ls (classic easy simple) or generic
p_controller0.set_type("ls")
p_controller0.set_nwfs([0])         # /!\
p_controller0.set_ndm([0,])       # /!\
p_controller0.set_maxcond(1500)     #     determines the nb of modes to be filtered
p_controller0.set_delay(Ts) # /!\ same delay in ms as in saxo.py
p_controller0.set_gain(0.3)


