import cvxpy as cp
import numpy as np

def rcone(x,y,z):
    v = 1/np.sqrt(2)
    V = np.array([[v,v],[v,-v]])
    T = np.block([[V,np.zeros((2,2))],[np.zeros((2,2)),np.eye(2)]])



if __name__ == '__main__':

    X_n = cp.Variable((2,1))
    Y_n = cp.vstack([1,cp.Variable(1)])
    XY_n = cp.vstack([X_n,Y_n])
    gam = cp.Variable((2,1))
    W1 = np.array([19,8])
    ZFy = np.array([[1.0000 + 0.0050j,1.0000 + 0.0000j],[1.0000 + 0.0100j,1.0000 + 0.0000j]])
    ZFx = np.array([[1.0000 + 0.0050j,1.0000 + 0.0000j],[1.0000 + 0.0100j,1.0000 + 0.0000j]])

    F_a = W1*ZFy@Y_n
    Pc = np.array([[1.2000 + 0.0050j],[1.2+0.01j]])
    G = np.array([[10000.0000 - 100.0050j],[10000-200.0050j]])

    PHI = 2 * cp.real(cp.multiply(cp.hstack([cp.multiply(G,ZFx), ZFy]), np.conj(Pc)))@XY_n - np.abs(Pc)**2


    v = 1/np.sqrt(2)
    V = np.array([[v,v],[v,-v]])
    T = np.block([[V,np.zeros((2,2))],[np.zeros((2,2)),np.eye(2)]])
    # T = np.insert(T,0,T[:,0],axis=1)
    # T = np.insert(T, 0, T[:, 0], axis=1)
    # T = np.insert(T, 0, T[:, 0], axis=1)
    plop = cp.vstack([PHI.T,gam.T,cp.real(F_a).T,cp.imag(F_a).T])
    Z = T@plop
    # soc_constraints = [cp.SOC(cp.real(F_a), cp.multiply(PHI,gam)),gam>=0,PHI>=0]
    # soc_constraints = [cp.SOC(PHI+gam,cp.vstack([2*F_a,PHI-gam]))]
    # soc_constraints = [cp.SOC(PHI[0] + gam[0], cp.vstack([2 * F_a[0], PHI[0] - gam[0]])),cp.SOC(PHI[1] + gam[1], cp.vstack([2 * F_a[1], PHI[1] - gam[1]]))]
    soc_constraints = [
        cp.SOC(PHI[i] + gam[i], cp.vstack([2 * F_a[i], PHI[i] - gam[i]])) for i in range(2)
    ]
    prob = cp.Problem(cp.Minimize(gam[0]+gam[1]), soc_constraints)

    prob.solve()