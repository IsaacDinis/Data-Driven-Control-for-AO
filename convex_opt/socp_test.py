import cvxpy as cp
import numpy as np

def rcone(x,y,z):
    v = 1/np.sqrt(2)
    V = np.array([[v,v],[v,-v]])
    T = np.block([[V,np.zeros((2,2))],[np.zeros((2,2)),np.eye(2)]])



if __name__ == '__main__':
    # n = 2
    # np.random.seed(2)
    # rand1 = np.random.randn(n)
    # rand2 = np.random.randn(n)
    # rand3 = np.random.randn(n)
    # rand4 = np.random.randn(1)
    # X = cp.Variable(n)
    # Y = cp.Variable(n)
    # ZFy = np.array([1+0.005j,1])
    # Yf = ZFy@Y
    # F_a = rand4*Yf
    # gamma = cp.Variable(1)
    # PHI = cp.hstack([cp.multiply(rand1,X)+rand3**2,cp.multiply(rand2,Y)+rand3**2])
    #
    #
    #
    # v = 1/np.sqrt(2)
    # V = np.array([[v,v],[v,-v]])
    # T = np.block([[V,np.zeros((2,2))],[np.zeros((2,2)),np.eye(2)]])
    # T = np.insert(T,0,T[:,0],axis=1)
    # T = np.insert(T, 0, T[:, 0], axis=1)
    # T = np.insert(T, 0, T[:, 0], axis=1)
    # plop = cp.hstack([PHI,gamma,cp.real(F_a),cp.imag(F_a)])
    # Z = T@plop
    # soc_constraints = [cp.SOC(Z[0], Z[2]),cp.SOC(Z[1], Z[3])]
    # prob = cp.Problem(cp.Minimize(gamma), soc_constraints)
    #
    # prob.solve()



    X_n = cp.Variable((2,1))
    Y_n = cp.vstack([1,cp.Variable(1)])
    XY_n = cp.vstack([X_n,Y_n])
    gam = cp.Variable(1)
    W1 = 19
    ZFy = np.array([1.0000 + 0.0050j,1.0000 + 0.0000j])
    ZFx = np.array([1.0000 + 0.0050j, 1.0000 + 0.0000j])

    F_a = W1*ZFy@Y_n
    Pc = 1.2000 + 0.0050j
    G = 10000.0000 - 100.0050j

    PHI = 2 * cp.real(cp.multiply(cp.hstack([cp.multiply(G,ZFx), ZFy]), np.conj(Pc)))@XY_n - np.abs(Pc)**2


    v = 1/np.sqrt(2)
    V = np.array([[v,v],[v,-v]])
    T = np.block([[V,np.zeros((2,2))],[np.zeros((2,2)),np.eye(2)]])
    # T = np.insert(T,0,T[:,0],axis=1)
    # T = np.insert(T, 0, T[:, 0], axis=1)
    # T = np.insert(T, 0, T[:, 0], axis=1)
    plop = cp.hstack([PHI,gam,cp.real(F_a),cp.imag(F_a)])
    Z = T@plop
    # soc_constraints = [cp.SOC(Z[2], Z[0])]
    soc_constraints = [Z[0]<=Z[2],Z[1]<=Z[3],gam>=0,Z[0]>=0,Z[1]>=0]
    prob = cp.Problem(cp.Minimize(gam), soc_constraints)

    prob.solve()