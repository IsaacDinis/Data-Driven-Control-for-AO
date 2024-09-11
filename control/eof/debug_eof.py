import astropy.io.fits as pfits
import numpy as np
from matplotlib import pyplot as plt

D = pfits.getdata('D.fits')
F = pfits.getdata('F.fits')
P = pfits.getdata('P.fits')
D = D[:,0:20000]
U, S, Vh = np.linalg.svd(D.T, full_matrices = False)
P = P[:,0:20000]
F_bis = P@U@(np.linalg.inv(np.diag(S))).T@Vh
# D /= np.linalg.norm(D,axis = 1).reshape(-1,1)
# P /= np.linalg.norm(P,axis = 1).reshape(-1,1)
# D_inv = np.linalg.pinv(D.T)
# F = (D_inv @ P.T).T
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