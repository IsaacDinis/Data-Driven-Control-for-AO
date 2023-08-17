#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import shesha.config as conf
import os
import atmosphere

simul_name = "2dm_1wfs"

Ts = 1./1000.
# loop
p_loop = conf.Param_loop()
p_loop.set_niter(5000)         # /?\ number of loops
p_loop.set_ittime(Ts) 

# Boston DMs                                # /?\ Choose the Boston DM
#boston_dm_file = 'Boston34x34_measured_IF.fits'
boston_dm_file = 'Boston34x34_flat.fits'                                                                                                                        
#boston_dm_file = 'Boston28x28_flat.fits'
#boston_dm_file = 'Boston24x24_flat.fits'
recommended_dm_pup = True                  # /?\ Choose using recommended DM pupil size

# geom
p_geom = conf.Param_geom()
p_geom.set_pupdiam(128)        # /!\ value to get WFS pixsize = 0.361 (SAXO truth = 0.36)
p_geom.set_zenithangle(0.)     # /!\ keep it 0 until we know what it does in COMPASS

# tel
p_tel = conf.Param_tel()
p_tel.set_diam(8.0)            # /!\  VLT diameter
# p_tel.set_type_ap('VLT-NoObs')       # /!\  VLT pupil
p_tel.set_cobs(0.14)  
p_tel.set_type_ap("VLT")       # VLT pupil
p_tel.set_spiders_type("four")
# p_tel.set_t_spiders(0.)
p_tel.set_t_spiders(0.00625)


p_atmos = conf.Param_atmos()
atm = atmosphere.TurbulenceProfile(single_layer=True)  # /?\ instantiate 1L or 35L atmosphere   
profileCondition = 'median'    # /?\ 'median', 'Q1', 'Q2', 'Q3' or 'Q4'
r0 = (atm.r0[atm.profile_conditions.index(profileCondition)]).value # /?\ take ESO r0 or..      
#r0 = 0.14                                                          # /?\ set your r0 value    
p_atmos.set_nscreens(len(atm.height.value))                         # /!\ set Nb of layers      
p_atmos.set_r0(r0)                                                  # /!\ set r0               
p_atmos.set_frac(atm.profile_data[profileCondition].tolist())       # /!\ set frac Cn2        
p_atmos.set_alt(atm.height.value)                                   # /!\ set heights          
p_atmos.set_windspeed(atm.wind_velocity.value)                      # /!\ set wind speeds      
p_atmos.set_winddir(atm.wind_dir.value)                             # /!\ set wind directions  
p_atmos.set_L0([25.] * len(atm.height.value))   


# target
p_target = conf.Param_target()
p_targets = [p_target]
p_target.set_xpos(0.)         # /!\ On axis
p_target.set_ypos(0.)         # /!\ On axis
p_target.set_Lambda(1.65)     # /!\ H Band
p_target.set_mag(1.)          # /!\

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
p_wfs0.set_gsmag(1.)
p_wfs0.set_optthroughput(1) # still unknown
p_wfs0.set_zerop(1e11)        # zero point for guide star magnitude
p_wfs0.set_noise(-1)         # EMCCD with < 0.1e- RON
p_wfs0.set_atmos_seen(1)      # /!\
p_wfs0.set_fstop("square")    # /!\
                              # Choose one spatial filter or none.
#p_wfs0.set_fssize(0.79412)   # 1.1*lambda/dSubap
#p_wfs0.set_fssize(0.8663)    # 1.2*lambda/dSubap
#p_wfs0.set_fssize(0.9385)    # 1.3*lambda/dSubap
#p_wfs0.set_fssize(1.0107)    # 1.4*lambda/dSubap
p_wfs0.set_fssize(1.0829)     # 1.5*lambda/dSubap
#p_wfs0.set_fssize(1.227275)  # 1.7*lambda/dSubap
#p_wfs0.set_fssize(1.44385)   # 2*lambda/dSubap

# dm
p_dm0 = conf.Param_dm()       # /!\
p_dm1 = conf.Param_dm()       # /!\
p_dms = [p_dm0, p_dm1]        # /!\
# p_dms = [p_dm0]        # /!\

p_dm0.set_type("pzt")         # /!\
p_dm0.set_alt(0.)             # /!\
p_dm0.set_unitpervolt(1.)     # /!\
p_dm0.set_thresh(-0.5)
p_dm0.set_file_influ_fits("HODM_gauss_fitSPARTA.fits") # /!\ to use a custom SAXO HO DM


p_dm1.set_type("pzt")           # /!\
p_dm1.set_thresh(-200)          # /!\ to get all Boston actuators
p_dm1.set_alt(0.)               # /!\
p_dm1.set_unitpervolt(1.)       # /!\
p_dm1.set_push4imat(1.0e-3)
p_dm1.set_file_influ_fits(boston_dm_file) # /!\ choice made at the begin. of this file
boston_dms = ['Boston24x24_flat.fits',
              'Boston28x28_flat.fits',
              'Boston34x34_flat.fits',
              'Boston34x34_measured_IF.fits']
if (boston_dm_file in boston_dms) and recommended_dm_pup:
    # recommanded pupil diam on DMs
    boston_rec_diam = [(24-3)*450e-6, (28-3)*450e-6,
                       (34-3)*400e-6, (34-3)*400e-6]
    p_dm1.set_diam_dm(boston_rec_diam[boston_dms.index(boston_dm_file)])
else:
    print('No recommended pupil size used on Boston DM')


# centroiders
p_centroider0 = conf.Param_centroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)     # /!\
p_centroider0.set_type("wcog")
p_centroider0.set_width(2) # size of the diffraction limit spot
p_centroider0.set_thresh(0)
# p_centroider0.set_type_fct("model")

# controllers
p_controller0 = conf.Param_controller()
p_controllers = [p_controller0]

p_controller0.set_type("generic")
p_controller0.set_nwfs([0])         # /!\
p_controller0.set_ndm([0,1])       # /!\
p_controller0.set_delay(1) # /!\ same delay in ms as in saxo.py
p_controller0.set_gain(0.3)

# coronagraph
# p_corono0 = conf.Param_corono()
# p_coronos = [p_corono0]

# p_corono0.set_type('SPHERE_APLC')
# p_corono0.set_dim_image(200)       # size of the science image in pixel
# p_corono0.set_wavelength_0(1.667)  # coronagraph wavelength in microns

# for a polychromatic coronagraph (optional), set these two parameters:
# p_corono0.set_nb_wav(3)           # number of simulated wavelength
# p_corono0.set_delta_wav(0.054)    # spectral bandwidth in micron
                                    # 0.054 μm = bandwidth of H3 IRDIS filter
                                    # this value is just an example, not a fixed paramete
