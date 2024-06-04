import cvxpy as cp
import numpy as np

def rcone(x,y,z):
    v = 1/np.sqrt(2)
    V = np.array([[v,v],[v,-v]])
    T = np.block([[V,np.zeros((2,2))],[np.zeros((2,2)),np.eye(2)]])



if __name__ == '__main__':
    n = 2
    np.random.seed(2)
    rand1 = np.random.randn(n)
    rand2 = np.random.randn(n)
    rand3 = np.random.randn(n)
    rand4 = np.random.randn(1)
    X = cp.Variable(n)
    Y = cp.Variable(n)
    ZFy = np.array([1+0.005j,1])
    Yf = ZFy@Y
    F_a = rand4*Yf
    gamma = cp.Variable(1)
    PHI = np.array([cp.multiply(rand1,X)+rand3**2,cp.multiply(rand2,Y)+rand3**2])



    v = 1/np.sqrt(2)
    V = np.array([[v,v],[v,-v]])
    T = np.block([[V,np.zeros((2,2))],[np.zeros((2,2)),np.eye(2)]])
    plop = cp.hstack([PHI,gamma,cp.real(F_a),cp.imag(F_a)])
    Z = T@plop