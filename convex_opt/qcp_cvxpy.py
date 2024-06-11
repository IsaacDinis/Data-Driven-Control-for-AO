import cvxpy as cp
import numpy as np

x = cp.Variable(1)
y = cp.Variable(1)
z = cp.Variable(1)

c1 = x + y + z == 1
c2 = x + y - z <= 0
c3 = cp.SOC(z+y,cp.vstack([2*x, z-y]))

constraints = [c1,c2,c3]
prob = cp.Problem(cp.Maximize(x), constraints)

prob.solve()