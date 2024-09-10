#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Friday 16th of September 2022

@author: Clementine Bechet for the SPHERE+ simulation group
"""

# /!\ : This symbol marks the parameters that must not be changed.
# /?\ : This symbol marks the parameters that must be chosen among a list.

import shesha.config as conf
import os
# import atmosphere
import numpy as np
simul_name = "1dm_1wfs"

Ts = 1./1000.
# loop
p_loop = conf.Param_loop()
p_loop.set_niter(5000)         # /?\ number of loops
p_loop.set_ittime(Ts) 


# geom
p_geom = conf.Param_geom()
p_geom.set_pupdiam(150)        # /!\ value to get WFS pixsize = 0.361 (SAXO truth = 0.36)
p_geom.set_zenithangle(0.)     # /!\ keep it 0 until we know what it does in COMPASS

# tel
p_tel = conf.Param_tel()
p_tel.set_diam(8.12)            # /!\  VLT diameter
# p_tel.set_type_ap('VLT-NoObs')       # /!\  VLT pupil
p_tel.set_cobs(0.16)  
p_tel.set_type_ap("VLT")       # VLT pupil
p_tel.set_spiders_type("four")
# p_tel.set_t_spiders(0.)
p_tel.set_t_spiders(0.00625)
# p_tel.set_t_spiders(0.0)


# atmos
p_atmos = conf.Param_atmos()
# p_atmos.set_r0(0.137)       # Fried parameters @ 500 nm
p_atmos.set_r0(500e-9/(4.85e-6*0.5))       # Fried parameters @ 500 nm
p_atmos.set_nscreens(1)    # Number of layers
p_atmos.set_frac([1.0])    # Fraction of atmosphere (100% = 1)
p_atmos.set_alt([0.0])     # Altitude(s) in meters
p_atmos.set_windspeed([6]) # wind speed of layer(s) in m/s
p_atmos.set_winddir([45])  # wind direction in degrees
p_atmos.set_L0([22])       # in meters
p_atmos.set_seeds(0)

# target
p_target = conf.Param_target()
p_targets = [p_target]
p_target.set_xpos(0.)         # /!\ On axis
p_target.set_ypos(0.)         # /!\ On axis
p_target.set_Lambda(0.75)     # /!\ H Band
p_target.set_mag(1.)          # /!\

# wfs
p_wfs0 = conf.Param_wfs(roket=False)
p_wfss = [p_wfs0]

p_wfs0.set_type("pyrhr")        # /!\ pyramid
p_wfs0.set_nxsub(60)            #     number of pixels
p_wfs0.set_fracsub(0.5)       #     threshold on illumination fraction for valid pixel
p_wfs0.set_Lambda(1.4)          #     wavelength
p_wfs0.set_gsmag(6.)
p_wfs0.set_zerop(1.e11)
p_wfs0.set_optthroughput(1)
p_wfs0.set_noise(-1)           #     readout noise
p_wfs0.set_xpos(0.)             # /!\ On axis
p_wfs0.set_ypos(0.)             # /!\ On axis

p_wfs0.set_pyr_ampl(0)
p_wfs0.set_pyr_npts(1) 
p_wfs0.set_pyr_pup_sep(p_wfs0.nxsub) # separation between the 4 images of the pyramid 
p_wfs0.set_fstop("round")
# p_wfs0.set_fssize(1.5)          # Size of the field stop
p_wfs0.set_fssize(3)  
p_wfs0.set_atmos_seen(1)        # /!\


p_centroider0 = conf.Param_centroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)           # /!\
p_centroider0.set_type("maskedpix")

# dm
p_dm0 = conf.Param_dm()       # /!\
p_dms = [p_dm0]        # /!\
# p_dms = [p_dm0, p_dm1]        # /!\



p_dm0.set_type("pzt")         # /!\
p_dm0.set_nact(48)
# p_dm0.set_thresh(0.485)        # /!\ to get the SAXO 1377 active actuators
p_dm0.set_thresh(0.)  
p_dm0.set_alt(0.)             # /!\
p_dm0.set_unitpervolt(1.)     # /!\
p_dm0.set_push4imat(0.180)    #     to displace ~ half a pixel
p_dm0.set_coupling(0.13)
# dead act = 305




# controllers
p_controller0 = conf.Param_controller()
p_controllers = [p_controller0]

p_controller0.set_type("generic")
p_controller0.set_nwfs([0])         # /!\
p_controller0.set_ndm([0])       # /!\
# p_controller0.set_ndm([0,1])       # /!\
p_controller0.set_delay(0) # /!\ same delay in ms as in saxo.py
p_controller0.set_gain(0.5)

p_corono0 = conf.Param_corono()
p_coronos = [p_corono0]

p_corono0.set_type('perfect')
p_corono0.set_dim_image(400)       # size of the science image in pixel
p_corono0.set_wavelength_0(0.75)  # coronagraph wavelength in microns
p_corono0.set_image_sampling(8)  # size of lambda/D in pixel