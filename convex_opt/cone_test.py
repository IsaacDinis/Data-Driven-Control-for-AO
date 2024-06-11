import cvxpy as cp
import numpy as np

x = cp.Variable(4)
z = cp.Variable(4)
y = cp.Variable(1)
c = np.array([-1.42,-0.15,-0.26,2.23])
a = np.array([-2.43,0.11,0.37,1.35])

soc_constraints = [cp.SOC(y,cp.multiply(z,x) + c)]
prob = cp.Problem(cp.Minimize(y),soc_constraints + [y >= 0])

prob.solve()
print(np.linalg.norm(np.multiply(a,x.value)+c))
# print(np.linalg.norm(a.T@x.value+d))
print(y.value)
