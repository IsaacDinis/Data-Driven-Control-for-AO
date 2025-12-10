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
import numpy as np
simul_name = "2dm_1wfs"

fs = 4000 #Hz
Ts = 1./fs

# loop
p_loop = conf.ParamLoop()
p_loop.set_niter(5000)         # /?\ number of loops
p_loop.set_ittime(Ts) 


# geom
p_geom = conf.ParamGeom()
p_geom.set_pupdiam(240)        # /!\ value to get WFS pixsize = 0.361 (SAXO truth = 0.36)
p_geom.set_zenithangle(0.)     # /!\ keep it 0 until we know what it does in COMPASS

# tel
p_tel = conf.ParamTel()
p_tel.set_diam(8.12)            # /!\  VLT diameter
# p_tel.set_type_ap('VLT-NoObs')       # /!\  VLT pupil
p_tel.set_cobs(0.16)  
p_tel.set_type_ap("VLT")       # VLT pupil
p_tel.set_spiders_type("four")
# p_tel.set_t_spiders(0.)
p_tel.set_t_spiders(0.00625)
# p_tel.set_t_spiders(0.0)


# atmos
p_atmos = conf.ParamAtmos()
# p_atmos.set_r0(0.137)       # Fried parameters @ 500 nm 0.75 arcs
# p_atmos.set_r0(0.2)       # Fried parameters @ 500 nm 0.75 arcs
# p_atmos.set_r0(500e-9/(4.85e-6*1.3))       # Fried parameters @ 500 nm
#p_atmos.set_r0(500e-9/(4.85e-6*0.88))       # Fried parameters @ 500 nm
#p_atmos.set_nscreens(1)    # Number of layers
#p_atmos.set_frac([1.0])    # Fraction of atmosphere (100% = 1)
#p_atmos.set_alt([0.0])     # Altitude(s) in meters
#p_atmos.set_windspeed([9.5]) # wind speed of layer(s) in m/s
#p_atmos.set_winddir([45])  # wind direction in degrees
#p_atmos.set_L0([50])       # in meters

#            nlayer         = 9; % # of layers
#            L0             = 25 ; % 25 [Nelly] % Turbulence external scale (m) % 22 = Eris sims
#            altitude       = [0.042 0.140 0.281 0.562 1.125 2.25 4.5 9  18]*1000; % Altitude of each layer in the atmosphere(m)
#            windSpeed      = [15    13    13    9     9     15   2    40   21]  ; % Wind speed on each layers in pure frozen flow
#            fractionnalR0  = [53.28 1.45  3.5   9.57  10.83 4.37 6.58 3.71 6.71]/100 ; % Cn2 profile coefficient for each layer
#            windDirection  = [38 34 54 42 57 48 -102 -83 -77]*pi/180; %2*pi*rand(1,param.nlayer) ; % Wind direction on each layer, compare to the ground
p_atmos.set_nscreens(9)    # Number of layers
p_atmos.set_r0(500e-9/(4.85e-6*0.88))       # Fried parameters @ 500 nm
p_atmos.set_frac([1.0])    # Fraction of atmosphere (100% = 1)
p_atmos.set_alt([0.042*1000, 0.140*1000, 0.281*1000, 0.562*1000, 1.125*1000, 2.25*1000, 4.5*1000, 9*1000,  18*1000])     # Altitude(s) in meters
p_atmos.set_windspeed([15,    13,    13,    9,     9,     15,   2,    40,   21]) # wind speed of layer(s) in m/s
p_atmos.set_winddir([38, 34, 54, 42, 57, 48, -102, -83, -77])  # wind direction in degrees
p_atmos.set_L0([50])       # in meters

p_atmos.set_seeds(123)

# target
p_target = conf.ParamTarget()
p_targets = [p_target]
p_target.set_xpos(0.)         # /!\ On axis
p_target.set_ypos(0.)         # /!\ On axis
p_target.set_Lambda(0.73)     # /!\ H Band
p_target.set_mag(1.)          # /!\

# wfs
p_wfs0 = conf.ParamWfs(roket=False)
p_wfss = [p_wfs0]

p_wfs0.set_type("pyrhr")        # /!\ pyramid
p_wfs0.set_nxsub(60)            #     number of pixels
p_wfs0.set_fracsub(0.5)       #     threshold on illumination fraction for valid pixel
p_wfs0.set_Lambda(1.4)          #     wavelength
photFlux = 4.31e7 # photon/m2/s
zerop = 1 # not used
gmag =  -2.5*np.log10(photFlux /zerop)
p_wfs0.set_gsmag(gmag)
p_wfs0.set_zerop(zerop)
p_wfs0.set_optthroughput(0.2)
p_wfs0.set_noise(0.5)           #     readout noise
p_wfs0.set_xpos(0.)             # /!\ On axis
p_wfs0.set_ypos(0.)             # /!\ On axis

p_wfs0.set_pyr_ampl(0)
p_wfs0.set_pyr_npts(1) 
p_wfs0.set_pyr_pup_sep(p_wfs0.nxsub) # separation between the 4 images of the pyramid 
p_wfs0.set_fstop("round")
p_wfs0.set_fssize(3)  
p_wfs0.set_atmos_seen(1)        # /!\

p_centroider0 = conf.ParamCentroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)           # /!\
p_centroider0.set_type("maskedpix")



# dm
p_dm0 = conf.ParamDm()       # /!\
p_dm1 = conf.ParamDm() 
p_dm2 = conf.ParamDm() 
     # /!\
p_dms = [p_dm0,p_dm1,p_dm2]        # /!\

p_dm0.set_type("pzt")         # /!\
nact = 41
p_dm0.set_nact(nact)
p_dm0.set_thresh(0.5)        # /!\ 
p_dm0.set_coupling(0.13)
p_dm0.set_alt(0.)             # /!\
p_dm0.set_unitpervolt(1.)     # /!\

p_dm1.set_type("pzt")         # /!\
p_dm1.set_thresh(-0.5)        # /!\ to get the SAXO 1377 active actuators
p_dm1.set_alt(0.)             # /!\
p_dm1.set_unitpervolt(1.)     # /!\
p_dm1.set_push4imat(0.180)    # to displace about half a pixel
p_dm1.set_file_influ_fits('HODM_gauss_fitSPARTA.fits') # /!\ to use a custom SAXO HO DM

# tip-tilt
p_dm2.set_type("tt")         # /!\
p_dm2.set_alt(0.)            # /!\
p_dm2.set_unitpervolt(1.)    # /!\
p_dm2.set_push4imat(0.18)    # to displace about half a pixel

# p_dm1.set_type("pzt")         # /!\
# # p_dm0.set_thresh(-0.1)        # /!\ to get the SAXO 1377 active actuators
# p_dm1.set_thresh(0.5)        # /!\ to get the SAXO 1377 active actuators
# p_dm1.set_alt(0.)             # /!\
# p_dm1.set_unitpervolt(1.)     # /!\
# p_dm1.set_push4imat(0.180)    #     to displace ~ half a pixel
# p_dm1.set_file_influ_fits("ristretto_43.fits") # /!\ to use a custom SAXO HO DM
# # p_dm1.set_file_influ_fits("HODM_gauss_fitSPARTA.fits") # /!\ to use a custom SAXO HO DM

# controllers
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

p_controller0.set_type("generic")
p_controller0.set_nwfs([0])         # /!\
p_controller0.set_ndm([0,1,2])       # /!\
if fs == 4000:
     p_controller0.set_delay(1) # /!\ 
else:
     p_controller0.set_delay(0.5)
p_controller0.set_gain(0)




#Sampling=240pixels

#Model atmosphere:

#            nlayer         = 9; % # of layers
#            L0             = 25 ; % 25 [Nelly] % Turbulence external scale (m) % 22 = Eris sims
#            altitude       = [0.042 0.140 0.281 0.562 1.125 2.25 4.5 9  18]*1000; % Altitude of each layer in the atmosphere(m)
#            windSpeed      = [15    13    13    9     9     15   2    40   21]  ; % Wind speed on each layers in pure frozen flow
#            fractionnalR0  = [53.28 1.45  3.5   9.57  10.83 4.37 6.58 3.71 6.71]/100 ; % Cn2 profile coefficient for each layer
#            windDirection  = [38 34 54 42 57 48 -102 -83 -77]*pi/180; %2*pi*rand(1,param.nlayer) ; % Wind direction on each layer, compare to the ground

#W#ave science: 730nm

#Threshold valid pixel WFS: 0.5

#Photons science: infini

#Nombre de modes controlles: 1190

#Couplage meca actuateurs: 13% (BMC)

#Retard: 2 frame

#Controlleur: OMGI vs DD4AO
