"""
atmosphere - Atmospheric parameters and functionality to be used in analysis code for SAXO+.
Courtesy of METIS team

author = 'Thomas Bertram' for METIS team
modified by C. Bechet 2023.03.16 for SAXO+ team.
"""

import numpy as np
import astropy.units as u
import pdb;

class TurbulenceProfile(object):
    """
    Atmospheric turbulence parameters for E-ELT AO analysis and simulations

    """

    def __init__(self, single_layer : bool = False):
        """
        Init of turbulence profile object

        Args:
            single_layer : boolean : If True, only one layer. Otherwise the 35-layer profile of
                                    ESO is used.

        """
        if single_layer is True:
            self.profile_data = np.array(
                [(1, 0, 10., 45., 1., 1., 1., 1., 1.)],
                dtype = [('layer', '<i4'), ('height', '<i4'), ('v_wind', '<f4'),
	                 ('d_wind', '<f4'), ('median', '<d'), ('Q1', '<d'), ('Q2', '<d'),
                         ('Q3', '<d'),('Q4', '<d')])
        else:
            
            # set a wind direction :: Cantalloube et al. ?High contrast imaging with ELT/METIS:
            # The wind driven halo, from SPHERE to METIS?. arXiv:1911.11241 (Nov. 2019),
            self.profile_data = np.array(
                [(1, 30, 5.5, 69.3952, 0.242, 0.226, 0.251, 0.255, 0.236),
                 (2, 90, 5.5, 69.72574, 0.12, 0.112, 0.116, 0.119, 0.131),
                 (3, 150, 5.1, 70.056274, 0.0968, 0.101, 0.0957, 0.0932, 0.0981),
                 (4, 200, 5.5, 70.33172, 0.059, 0.064, 0.0584, 0.0557, 0.0577),
                 (5, 245, 5.6, 70.57961, 0.0473, 0.0415, 0.037, 0.045, 0.0658),
                 (6, 300, 5.7, 70.8826, 0.0473, 0.0415, 0.037, 0.045, 0.0658),
                 (7, 390, 5.8, 71.3784, 0.0473, 0.0415, 0.037, 0.045, 0.0658),
                 (8, 600, 6., 72.6508, 0.0473, 0.0415, 0.037, 0.045, 0.0658),
                 (9, 1130, 6.5, 76.79038, 0.0399, 0.031, 0.0325, 0.0419, 0.0540),
                 (10, 1880, 7., 93.52274, 0.0324, 0.0226, 0.0347, 0.0404, 0.0320),
                 (11, 2630, 7.5, 185.61185, 0.0162, 0.0113, 0.0174, 0.0202, 0.0160),
                 (12, 3500, 8.5, 267.1189, 0.026, 0.0221, 0.03, 0.0304, 0.0218),
                 (13, 4500, 9.5, 281.4172, 0.0156, 0.0133, 0.018, 0.0182, 0.0131),
                 (14, 5500, 11.5, 284.26138, 0.0104, 0.0088, 0.012, 0.0121, 0.0087),
                 (15, 6500, 17.5, 286.76254, 0.01, 0.0147, 0.013, 0.0086, 0.0037),
                 (16, 7500, 23., 288.81754, 0.012, 0.0177, 0.0156, 0.0103, 0.0045),
                 (17, 8500, 26., 281.50632, 0.004, 0.0059, 0.0052, 0.0034, 0.0015),
                 (18, 9500, 29., 274.19513, 0.014, 0.0206, 0.0182, 0.012, 0.0052),
                 (19, 10500, 32., 276.93207, 0.013, 0.0192, 0.017, 0.0111, 0.0049),
                 (20, 11500, 27., 277.9622, 0.007, 0.0103, 0.0091, 0.006, 0.0026),
                 (21, 12500, 22., 277.63654, 0.016, 0.023, 0.0187, 0.0143, 0.0080),
                 (22, 13500, 14.5, 275.8745, 0.0259, 0.0375, 0.0303, 0.0231, 0.0129),
                 (23, 14500, 9.5, 276.17603, 0.019, 0.0276, 0.0223, 0.017, 0.0095),
                 (24, 15500, 6.3, 278.88297, 0.0099, 0.0143, 0.0115, 0.0088, 0.0049),
                 (25, 16500, 5.5, 281.443, 0.0062, 0.0089, 0.0072, 0.0055, 0.0031),
                 (26, 17500, 6., 275.80304, 0.004, 0.0058, 0.0047, 0.0036, 0.0020),
                 (27, 18500, 6.5, 270.16306, 0.0025, 0.0036, 0.003, 0.0022, 0.0012),
                 (28, 19500, 7., 272.53885, 0.0022, 0.0031, 0.0025, 0.0019, 0.0010),
                 (29, 20500, 7.5, 276.1623, 0.0019, 0.0027, 0.0022, 0.0017, 0.0009),
                 (30, 21500, 8., 275.62222, 0.0014, 0.002, 0.0016, 0.0012, 0.0007),
                 (31, 22500, 8.5, 274.08298, 0.0011, 0.0016, 0.0013, 0.001, 0.0006),
                 (32, 23500, 9., 272.54376, 0.0006, 0.0009, 0.0007, 0.0006, 0.0003),
                 (33, 24500, 9.5, 262.65988, 0.0009, 0.0012, 0.0011, 0.0008, 0.0005),
                 (34, 25500, 10., 247.50516, 0.0005, 0.0007, 0.0006, 0.0004, 0.0002),
                 (35, 26500, 10., 232.35045, 0.0004, 0.0006, 0.0005, 0.0004, 0.0002)],
                dtype=[('layer', '<i4'), ('height', '<i4'), ('v_wind', '<f4'),
                       ('d_wind', '<f4'), ('median', '<d'), ('Q1', '<d'), ('Q2', '<d'),
                       ('Q3', '<d'), ('Q4', '<d')])

        self.wind_velocity = self.profile_data['v_wind'] * u.m / u.s

        self.wind_dir = self.profile_data['d_wind'] * u.deg 
            
        self.height = self.profile_data['height'] * u.m
            
        self.profile_conditions = ['median', 'Q1', 'Q2', 'Q3', 'Q4']
        
        self.fractional = self.profile_data[self.profile_conditions]
        
        self.reference_wl = 0.5 * u.micron
        
        self.r0 = np.array((0.157, 0.234, 0.178, 0.139, 0.097)) * u.m
        
        self.alpha = np.array((1.052, 0.925, 0.968, 1.079, 1.370),
                              dtype=[('median', '<f4'), ('Q1', '<f4'), ('Q2', '<f4'),
                                     ('Q3', '<f4'), ('Q4', '<f4')])

            
    @u.quantity_input
    def coherence_length(self, wavelength: u.micron = None,
                         zenith_angle: u.deg = 0. *u.deg,
                         profile_condition = None):
        """Returns the turbulence coherence length (Fried parameter) r_0
        
        The turbulence coherence length, as provided in ESO-258292 is 
        available for the median profile condition as well as for each quartile 
        of the profile condition determined in long term cn^2 measurements for 
        Armazones. It can be scaled to the wavefront sensing wavelength and the
        zenith angle. 
         
        If no wavelength is specified, the reference wavelength (0.5 micron) is assumed. 
        
        If no zenith angle is specified, zenith pointing is assumed.
        
        Args:
            wavelength (length unit): (optional) wavefront sensing wavelength. 
            zenithAngle (angle unit): (optional) zenith angle. 
            profileCondition (str): (optional) 'median', 'Q1', 'Q2', 'Q3' or 'Q4'

        Returns:
            ndarray with turbulence coherence lengths for all profile conditions or
            quantity with turbulence coherence length in m for the specified 
            profile condition
        """
        if profile_condition is None:
            temp_r0 = self.r0
        else:
            if profile_condition not in self.profile_conditions:
                raise (ProfileConditionSpecifierError)
            temp_r0 = self.r0[self.profile_conditions.index(profile_condition)]

        if wavelength is None:
            wavelength = self.reference_wl
            
        return temp_r0 * (wavelength.to(u.micron) / self.reference_wl) ** (6. / 5.) * (np.cos(zenith_angle) ** (3. / 5.))


    @u.quantity_input
    def coherence_time(self, wavelength: u.micron = None,
                       zenith_angle: u.deg = 0. *u.deg,
                       profile_condition = None):
        """Returns the turbulence coherence time tau_0 = 0.31 r0 / V.                        

        If no wavelength is specified, the reference wavelength (0.5 micron) is assumed.        

        If no zenith angle is specified, zenith pointing is assumed.                             
        Args:                                                                                    
            wavelength (length unit): (optional) wavefront sensing wavelength.                   
            zenithAngle (angle unit): (optional) zenith angle.                                   
            profileCondition (str): 'median' (default), 'Q1', 'Q2', 'Q3' or 'Q4'                
                                                                                                        Returns:                                                                                 
            ndarray with turbulence coherence lengths for all profile conditions or              
            quantity with turbulence coherence length in m for the specified                     
            profile condition                                                                    
        """
        if profile_condition is None:
            profile_condition = 'median'

        if profile_condition not in self.profile_conditions:
            raise (ProfileConditionSpecifierError)

        r0 = self.r0[self.profile_conditions.index(profile_condition)]

        if wavelength is None:
            wavelength = self.reference_wl

        wind_v = self.wind_velocity / np.cos(zenith_angle) # rescale V w.r.t. zenith angle
        
        v = np.sum(wind_v**(5./3) * self.fractional[profile_condition])**(3./5)
    
        return 0.31 * r0 * (wavelength.to(u.micron) / self.reference_wl) ** (6. / 5.) * (np.cos(zenith_angle) ** (3. / 5.)) / v   


    @u.quantity_input
    def heights(self, zenith_angle: u.deg = 0. * u.deg):
        """Returns the layer heights along the line of sight 
        The layer altitudes provided by ESO can be scaled along the line of sight with the 
        zenith angle.
        If no zenith angle is specified, zenith pointing is assumed.                          

        Args:
            zenithAngle (angle unit): (optional) zenith angle.                                 
                Defaults to 0 deg.                                                             

        Returns:                                                                               
          ndarray with heights in m

        """
        return self.profile_data['height'].tolist() / np.cos(zenith_angle)


class Error(Exception):
    """Base class for exceptions in the saxoplus atmosphere module."""
    pass


class ProfileConditionSpecifierError(Error):
    def __str__(self):
        return "Specification of the profile condition is wrong. Use 'median', 'Q1', 'Q2' ,'Q3', or 'Q4'"
