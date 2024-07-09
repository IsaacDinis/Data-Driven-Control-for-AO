import cvxpy as cp
import numpy as np

x = cp.Variable(1)
y = cp.Variable(1)
z = cp.Variable(1)

c1 = x + y + z == 1
c2 = x + y - z <= 0
c3 = cp.SOC(z+y,cp.vstack([2*x, z-y]))

v = 1 / np.sqrt(2)
V = np.array([[v, v], [v, -v]])
T = np.block([[V, np.zeros((2, 2))], [np.zeros((2, 2)), np.eye(2)]])
plop = cp.vstack([y,z,np.sqrt(2)*x,0])
Z = T@plop
c3 = cp.SOC(Z[0],Z[1:])
constraints = [c1,c2,c3]
prob = cp.Problem(cp.Maximize(x), constraints)

prob.solve()