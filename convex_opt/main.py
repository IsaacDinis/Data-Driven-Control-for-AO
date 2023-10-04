# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import cvxpy as cp
import numpy as np
def H2_opt_1_freq():
    n = 2
    G = 0.3090 - 0.9511j
    X_c = np.array([0.5,0])
    Y_c = np.array([1,-1])
    # X_c = 1.0000e+03 + 3.1416e+05j
    # Y_c = -1.0000e+03 + 6.2832e+05j
    X = cp.Variable((1,n))
    Y = cp.Variable((1,n))
    # Z = np.array([0.809016994374948 + 0.587785252292473j,1])
    # X_c *= Z
    # Y_c *= Z
    W1 = 0.6
    P = Y+G*X
    Pc = Y_c+G*X_c
    gamma_2 = cp.Variable((1,1))
    # TODO augment dimension variable en dessous
    dummy = cp.conj(P) @ Pc #+ cp.conj(Pc) @ P.T - cp.conj(Pc) @ Pc
    # cons = np.array([[gamma_2,W1*Y],[cp.conj(0.6*Y),cp.conj(P)*Pc+cp.conj(Pc)*P-cp.conj(Pc)*Pc]])
    cons = cp.vstack([cp.hstack([gamma_2, W1*Y]), cp.hstack([cp.conj(W1*Y), cp.conj(P)@Pc+cp.conj(Pc)@P.T-cp.conj(Pc)@Pc])])
    constraints = [cons >= 0 , gamma_2 >= 0]
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
