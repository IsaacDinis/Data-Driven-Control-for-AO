import numpy as np
h = np.array([1.,1.])
# F = np.array([2.865018208900050,-1.86505369404750])
F = np.array([1.10403652271300,-0.104155503015019])
for i in range(1000000):
    command = F@h
    h = np.roll(h,1)
    h[0] = command
print(h)
