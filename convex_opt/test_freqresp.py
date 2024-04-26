import control as ct
import numpy as np
Ts = 1/3000
szy = 2
szx = 2
W = 15
z = ct.tf([1,0],[0,1], dt = Ts)
z_ = ct.freqresp(z, W)
z_ = ct.freqresp(z, W).fresp.squeeze()
Zy = z_**np.arange(szy-1,-1,-1)
Zx = z_**np.arange(szy-1,-1,-1)