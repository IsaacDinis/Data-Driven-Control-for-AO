import numpy as np
import astropy.io.fits as pfits
from matplotlib import pyplot as plt

slopes_array = pfits.getdata('slopes_array.fits').squeeze()
command_array = pfits.getdata('command_array.fits').squeeze()
phase_array = pfits.getdata('phase_array.fits').squeeze()
phase_array_corr = pfits.getdata('phase_array_corr.fits').squeeze()