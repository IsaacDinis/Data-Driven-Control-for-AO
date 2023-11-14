#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Friday 16th of September 2022

@author: Clementine Bechet for the SPHERE+ simulation group
"""

# /!\ : This symbol marks the parameters that must not be changed.
# /?\ : This symbol marks the parameters that must be chosen among a list.

import shesha.config as conf
import numpy as np
import os


simul_name = "saxo"

# /?\ Choose a first stage loop frequency. (you may create your own pair)
saxo_frequency_delay_pair = [1380.,1.15] # from Cantalloube et al.

# loop
p_loop = conf.Param_loop()
p_loop.set_niter(5000)                                # number of loops
p_loop.set_ittime(1./saxo_frequency_delay_pair[0])   # first stage frequency 
p_loop.set_devices([0, 1, 2, 3])

# geom
p_geom = conf.Param_geom()
p_geom.set_pupdiam(400)        # /!\ value to get WFS pixsize = 0.338 (SAXO truth = 0.36)
p_geom.set_zenithangle(0.)

# tel
p_tel = conf.Param_tel()
p_tel.set_diam(8.0)            # /!\  VLT diameter
p_tel.set_cobs(0.14)           # /!\  central obstruction
p_tel.set_type_ap("VLT")       # /!\  VLT pupil
p_tel.set_t_spiders(0.00625)   # /!\  spider width = 5 cm

# atmos
p_atmos = conf.Param_atmos()
p_atmos.set_r0(0.15)       # Fried parameters @ 500 nm
p_atmos.set_nscreens(1)    # Number of layers
p_atmos.set_frac([1.0])    # Fraction of atmosphere (100% = 1)
p_atmos.set_alt([0.0])     # Altitude(s) in meters
p_atmos.set_windspeed([8]) # wind speed of layer(s) in m/s
p_atmos.set_winddir([45])  # wind direction in degrees
p_atmos.set_L0([25])       # in meters

# target
p_target = conf.Param_target()
p_targets = [p_target]
p_target.set_xpos(0.)         # /!\ On axis
p_target.set_ypos(0.)         # /!\ On axis
p_target.set_Lambda(1.65)     # /!\ H Band
p_target.set_mag(6.)          # /!\

# wfs
p_wfs0 = conf.Param_wfs()
p_wfss = [p_wfs0]

p_wfs0.set_type("sh")         # /!\ Shack-Hartmann
p_wfs0.set_nxsub(40)          # /!\ nb of sub-apertures.
p_wfs0.set_npix(6)            # /!\ nb of pixels / sub-aperture.
p_wfs0.set_pixsize(0.36)      # /!\ Shannon at 700nm. No exact reference found
p_wfs0.set_fracsub(0.5)       # /!\ Select 1240 subapertures.
p_wfs0.set_xpos(0.)           # /!\ On axis
p_wfs0.set_ypos(0.)           # /!\ On axis
p_wfs0.set_Lambda(0.7)        # /!\ SAXO SH bandwidth : [475, 900] nm
p_wfs0.set_gsmag(6.)
p_wfs0.set_optthroughput(0.5) # still unknown
p_wfs0.set_zerop(1e11)        # zero point for guide star magnitude
p_wfs0.set_noise(0.1)         # EMCCD with < 0.1e- RON
p_wfs0.set_atmos_seen(1)      # /!\
p_wfs0.set_fstop("square")    # /!\
                              # Choose one spatial filter or none.
p_wfs0.set_fssize(1.0829)     # 1.5*lambda/dSubap

# dm
p_dm0 = conf.Param_dm()       # /!\
p_dm1 = conf.Param_dm()       # /!\
p_dms = [p_dm0, p_dm1]        # /!\

p_dm0.set_type("pzt")         # /!\
p_dm0.set_thresh(-0.5)        # /!\ to get the SAXO 1377 active actuators
p_dm0.set_alt(0.)             # /!\
p_dm0.set_unitpervolt(1.)     # /!\
p_dm0.set_push4imat(0.180)    #     to displace ~ half a pixel
p_dm0.set_file_influ_fits("HODM_gauss_fitSPARTA.fits") # /!\ to use a custom SAXO HO DM

# tip-tilt
p_dm1.set_type("tt")         # /!\
p_dm1.set_alt(0.)            # /!\
p_dm1.set_unitpervolt(1.)    # /!\
p_dm1.set_push4imat(0.18)    #     to displace about half a pixel

# centroiders
p_centroider0 = conf.Param_centroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)     # /!\
p_centroider0.set_type("wcog")
p_centroider0.set_width(2) # size of the diffraction limit spot
p_centroider0.set_thresh(0)

# controllers
p_controller0 = conf.Param_controller()
p_controllers = [p_controller0]

p_controller0.set_type("generic")   # /?\ ls (classic easy simple) or generic
#p_controller0.set_type("ls")
p_controller0.set_nwfs([0])         # /!\
p_controller0.set_ndm([0, 1])       # /!\
p_controller0.set_maxcond(1500)     #     determines the nb of modes to be filtered
p_controller0.set_delay(0) # /!\ see at the head of the file.
p_controller0.set_gain(0.3)
