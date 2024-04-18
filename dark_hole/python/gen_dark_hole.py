from hcipy import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import astropy.io.fits as pfits
def plot_vapp(vAPP, prop):
    '''Plot the phase pattern and PSF of a vAPP

    Parameters
    ---------
    vAPP : Wavefront
        The wavefront of a vAPP mask, containing the vAPP pattern as phase and
        the telescope pupil as amplitude
    prop : Function
        A propagator function that propagates the wavefront to a focal plane
    '''

    # Plotting the phase pattern and the PSF
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    im1 = imshow_field(vAPP.phase, mask=vAPP.amplitude, cmap='RdBu')

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax=cax, orientation='vertical')

    ax2 = fig.add_subplot(122)
    im2 = imshow_field(np.log10(prop(vAPP).intensity/np.max(prop(vAPP).intensity)),vmin = -6, cmap='inferno')

    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im2, cax=cax, orientation='vertical')
    plt.show()

# m_size = 404 #saxo
m_size = 416 #saxo+
s_size = 400
pad_size = int((m_size-s_size)/2)

pupil_grid = make_pupil_grid(s_size)
focal_grid = make_focal_grid(4, 20)

prop = FraunhoferPropagator(pupil_grid, focal_grid)

aperture = make_vlt_aperture(True, with_spiders= False)
# aperture = make_obstructed_circular_aperture(8.0, 0., num_spiders=0, spider_width=0.00625)

telescope_pupil = aperture(pupil_grid)

imshow_field(telescope_pupil, cmap='gray')
plt.show()

contrast_level = 1e-7
dark_zone = (make_circular_aperture(12)(focal_grid)).astype(bool)*(focal_grid.x>4)
contrast = focal_grid.ones()
contrast[dark_zone] = contrast_level
imshow_field(np.log10(contrast))
plt.colorbar()
plt.show()

# Setting up the vAPP calculation parameters.
num_iterations = 80
wavefront = Wavefront(telescope_pupil, 1)

# Generate the vAPP pattern.
vAPP = generate_app_keller(wavefront, prop, contrast, num_iterations, beta = 1)

plot_vapp(vAPP, prop)
dark_hole_pad = np.pad(vAPP.phase.shaped,pad_size)
dark_hole = np.dstack([dark_hole_pad])
pfits.writeto('dark_hole.fits', dark_hole, overwrite = True)
