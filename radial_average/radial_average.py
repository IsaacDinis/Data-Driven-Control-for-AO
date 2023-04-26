from astropy.io import fits
from matplotlib import pyplot as plt

fits_image_filename = 'dummy.fits'

hdul = fits.open(fits_image_filename)
image = hdul[0].data
plt.figure(1)
plt.imshow(image)
plt.show()