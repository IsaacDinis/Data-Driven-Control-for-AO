#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Friday 7th of February 2020
@author: nour
"""

import shesha.config as conf


# loop
p_loop = conf.Param_loop()
p_loop.set_niter(5000)          # number of loop iterations
p_loop.set_ittime(1./1000.)     # assuming loop at 1.38 kHz

# geom
p_geom = conf.Param_geom()
p_geom.set_pupdiam(400)
p_geom.set_zenithangle(0.)

# tel
p_tel = conf.Param_tel()
p_tel.set_diam(8.0)            # VLT diameter
p_tel.set_cobs(0.14)           # central obstruction
p_tel.set_type_ap("VLT")       # VLT pupil
p_tel.set_spiders_type("four")
p_tel.set_t_spiders(0.00625)

# atmos
p_atmos = conf.Param_atmos()
p_atmos.set_r0(0.15)       # Fried parameters @ 500 nm
p_atmos.set_nscreens(1)    # Number of layers
p_atmos.set_frac([1.0])    # Fraction of atmosphere (100% = 1)
p_atmos.set_alt([0.0])     # Altitude(s) in meters
p_atmos.set_windspeed([34]) # wind speed of layer(s) in m/s
p_atmos.set_winddir([0])  # wind direction in degrees
p_atmos.set_L0([25])       # in meters

# target
p_target = conf.Param_target()
p_targets = [p_target]

p_target.set_xpos(0.)
p_target.set_ypos(0.)
p_target.set_Lambda(1.65) # H Band
p_target.set_mag(1.)

# wfs
p_wfs0 = conf.Param_wfs()
p_wfss = [p_wfs0]

p_wfs0.set_type("sh")
p_wfs0.set_nxsub(40)
p_wfs0.set_npix(6)
p_wfs0.set_pixsize(0.3)
p_wfs0.set_fracsub(0.5)
p_wfs0.set_xpos(0.)
p_wfs0.set_ypos(0.)
p_wfs0.set_Lambda(0.7)        # SAXO SH bandwidth : [475, 900] nm
p_wfs0.set_gsmag(1.)
p_wfs0.set_optthroughput(0.1) # still unknown
p_wfs0.set_zerop(1e11)        # zero point for guide star magnitude
p_wfs0.set_noise(-1)         # -1 = No noise
p_wfs0.set_atmos_seen(1)
p_wfs0.set_fstop("round")
#p_wfs0.set_fssize(0.79412)   # 1.1*lambda/dSubap
#p_wfs0.set_fssize(0.8663)    # 1.2*lambda/dSubap
#p_wfs0.set_fssize(0.9385)    # 1.3*lambda/dSubap
#p_wfs0.set_fssize(1.0107)    # 1.4*lambda/dSubap
p_wfs0.set_fssize(1.0829)     # 1.5*lambda/dSubap
#p_wfs0.set_fssize(1.227275)  # 1.7*lambda/dSubap
#p_wfs0.set_fssize(1.44385)   # 2*lambda/dSubap

# dm (waiting for the custom HODM)
p_dm0 = conf.Param_dm()
p_dm1 = conf.Param_dm()
p_dms = [p_dm0, p_dm1]

p_dm0.set_type("pzt")
nact = 41
p_dm0.set_nact(nact)
p_dm0.set_alt(0.)
p_dm0.set_thresh(0.25)
p_dm0.set_coupling(0.3)
p_dm0.set_unitpervolt(1.)
p_dm0.set_push4imat(0.001)
#p_dm0.set_margin_out(0.3)

# tip-tilt
p_dm1.set_type("tt")
p_dm1.set_alt(0.)
p_dm1.set_unitpervolt(1.)
p_dm1.set_push4imat(0.005)

# centroiders
p_centroider0 = conf.Param_centroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)
p_centroider0.set_type("bpcog")
p_centroider0.set_nmax(6) # nb of brightest pixels
# p_centroider0.set_type("corr")
# p_centroider0.set_type_fct("model")

# controllers
p_controller0 = conf.Param_controller()
p_controllers = [p_controller0]

p_controller0.set_type("generic") # ls (classic easy simple) or generic
p_controller0.set_nwfs([0])
p_controller0.set_ndm([0, 1])
p_controller0.set_maxcond(880)   # what determines the number of modes to be filtered
p_controller0.set_delay(1)
p_controller0.set_gain(0.3)
#p_controller0.set_nstates(6)