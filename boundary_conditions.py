#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import sys
import os

if __name__ == "__main__":
    constants_arg = "orchard"
    args = sys.argv
    if len(args) >= 2:
        constants_arg = args[1]

    os.chdir('C++')
    p = subprocess.run(f"./boundary_test {constants_arg}", shell=True)
    
    os.chdir('..')
    A = np.loadtxt("boundary_test.txt", delimiter=',')
    x = A[:, 0]
    cu = A[:, 1]
    cv = A[:, 2]

    fig, ax = plt.subplots()


    plt.loglog(x, cu, 'ro-')
    plt.loglog(x, cv, 'bo-')
    plt.legend(["O2", "CO2"])
    plt.xlabel("Maximum triangle size (m^2)")
    plt.ylabel("Maximum absolute error in middle of intervals")
    plt.title("Maximum error in the middle of the intervals of the boundary")
    plt.show()
