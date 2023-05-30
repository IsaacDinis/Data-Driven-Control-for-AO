# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from cvxopt import matrix, solvers
import numpy

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    A = matrix([ [-1.0, -1.0, 0.0, 1.0], [1.0, -1.0, -1.0, -2.0] ])
    b = matrix([1.0, -2.0, 0.0, 4.0])
    c = matrix([2.0, 1.0])
    sol = solvers.lp(c, A, b)
    print(sol['x'])
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
