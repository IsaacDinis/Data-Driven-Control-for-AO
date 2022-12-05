import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

fs = 2000

w_DM = 300*2*np.pi
s_DM = 0.1
DM_c = signal.lti([w_DM**2], [1,2*s_DM*w_DM, w_DM**2])

filts = DM_c
filtz = signal.lti(*signal.bilinear(filts.num, filts.den, fs))
wz, hz = signal.freqz(filtz.num, filtz.den)
ws, hs = signal.freqs(filts.num, filts.den, worN=fs*wz)

plt.semilogx(wz*fs/(2*np.pi), 20*np.log10(np.abs(hz).clip(1e-15)),label=r'$|H_z(e^{j \omega})|$')
plt.semilogx(wz*fs/(2*np.pi), 20*np.log10(np.abs(hs).clip(1e-15)),  label=r'$|H(j \omega)|$')
plt.legend()
plt.xlabel('Frequency [Hz]')
plt.ylabel('Magnitude [dB]')
plt.grid(True)
plt.show()