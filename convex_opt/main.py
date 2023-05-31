# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import cvxpy as cp
import numpy as np

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    m = 20
    n = 15
    np.random.seed(1)
    A = np.random.randn(m, n)
    b = np.random.randn(m)

    # Define and solve the CVXPY problem.
    x = cp.Variable(n)
    cost = cp.sum_squares(A @ x - b)
    prob = cp.Problem(cp.Minimize(cost))
    prob.solve()

    # Print result.
    print("\nThe optimal value is", prob.value)
    print("The optimal x is")
    print(x.value)
    print("The norm of the residual is ", cp.norm(A @ x - b, p=2).value)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
