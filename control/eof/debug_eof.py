import astropy.io.fits as pfits
import numpy as np
from matplotlib import pyplot as plt

D = pfits.getdata('D.fits')
F = pfits.getdata('F.fits')
P = pfits.getdata('P.fits')
U, S, Vh = np.linalg.svd(D.T, full_matrices = False)
F_bis = P@U@(np.linalg.inv(np.diag(S))).T@Vh
A = F@D
plt.figure()
plt.plot(F[0,:])
plt.show()

plt.figure()
plt.imshow(F,interpolation = None)
plt.show()

plt.figure()
plt.plot(P[0,:])
plt.show()

plt.figure()
plt.imshow(P@P.T,interpolation = None)
plt.show()