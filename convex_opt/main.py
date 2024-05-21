# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import cvxpy as cp
import numpy as np
import control as ct
def H2_opt_1_freq():
    Fx = 1
    Fy = 1
    n = 2
    # G = 1.0000 - 0.0050j
    X_c = np.array([0.2,0]).reshape(2, 1)
    Y_c = np.array([1,-1]).reshape(2, 1)
    X = cp.Variable((n,1))
    Y = cp.Variable((n,1))
    X_n = X+X_c
    Y_n = Y+Y_c
    XY_n = cp.vstack((X_n, Y_n))
    Ts = 1 / 3000
    szy = 2
    szx = 2
    W = 15
    z = ct.tf([1, 0], [0, 1], dt=Ts)
    G = ct.freqresp(1/z,W).fresp.squeeze()
    z_ = ct.freqresp(z, W).fresp.squeeze()
    Zy = z_ ** np.arange(szy - 1, -1, -1)
    Zx = z_ ** np.arange(szy - 1, -1, -1)
    ZFx = Zx * Fx
    ZFy = Zy * Fy
    Xcs = Zx@(X_c)
    Ycs = Zy @ (Y_c)
    Xc = Xcs* Fx
    Yc = Ycs * Fy

    W1 = 0.6
    P = Y+G*X
    Pc = Yc+G*Xc
    Cp = np.concatenate([G*ZFx, ZFy],axis= 0)
    x1_1 = -cp.multiply(cp.conj(Pc),Pc)
    x1_2 = 2*cp.real(cp.multiply(Cp, cp.conj(Pc)))
    x1 = x1_2@XY_n + x1_1
    gamma_2 = cp.Variable((1,1))
    x2 = gamma_2
    x3_d = cp.multiply(W1,ZFy)
    x3 = cp.hstack([cp.real(x3_d)@Y_n,cp.imag(x3_d)@Y_n])
    cons = cp.hstack([x1,x1,x3])
    # TODO augment dimension variable en dessous
    # dummy = cp.real(cp.conj(P)@Pc)+cp.real(cp.conj(Pc)@P.T)-cp.real(cp.conj(Pc)@Pc)
    # plop = cp.reshape(dummy,(1,1))
    # cons = np.array([[gamma_2,W1*Y],[cp.conj(0.6*Y),cp.conj(P)*Pc+cp.conj(Pc)*P-cp.conj(Pc)*Pc]])
    # cons = cp.vstack([cp.hstack([gamma_2, W1*Y]), cp.hstack([cp.conj(W1*Y), plop])])
    constraints = [cons >> 0 , gamma_2 >> 0]
    objective = cp.Minimize(cp.abs(gamma_2))

    prob = cp.Problem(objective, constraints)
    print("Optimal value", prob.solve())

def H2_opt():
    X_c = np.array([0.5,0])
    Y_c = np.array([1,-1])
def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # m = 20
    # n = 15
    # np.random.seed(1)
    # A = np.random.randn(m, n)
    # b = np.random.randn(m)
    #
    # # Define and solve the CVXPY problem.
    # x = cp.Variable(n)
    # cost = cp.sum_squares(A @ x - b)
    # prob = cp.Problem(cp.Minimize(cost))
    # prob.solve()
    #
    # # Print result.
    # print("\nThe optimal value is", prob.value)
    # print("The optimal x is")
    # print(x.value)
    # print("The norm of the residual is ", cp.norm(A @ x - b, p=2).value)

    H2_opt_1_freq()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
