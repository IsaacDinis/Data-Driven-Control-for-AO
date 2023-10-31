import numpy as np
import astropy.io.fits as pfits
from matplotlib import pyplot as plt

bump = pfits.getdata('../compass/results/bump_5/corono.fits')
no_bump = pfits.getdata('../compass/results/no_bump_5/corono.fits')
plt.figure()
coroimgSampl = 1/8
nimg = np.shape(bump)[0]
coroX = np.arange(-nimg // 2 * coroimgSampl,
                  nimg // 2 * coroimgSampl, coroimgSampl)
im = plt.pcolormesh(coroX, coroX, no_bump-bump)
plt.colorbar(im)
plt.xlabel(r'x [$\lambda$ / D]')
plt.ylabel(r'y [$\lambda$ / D]')
plt.title('corono image')
plt.show()