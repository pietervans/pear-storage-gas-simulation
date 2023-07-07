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
    p = subprocess.run(f"./convergence_test {constants_arg}", shell=True)

    
    os.chdir('..')
    A = np.loadtxt("convergence_test.txt", delimiter=',')
    sizes = A[:, 0]
    cu = A[:, 1]
    cv = A[:, 2]

    fig, ax = plt.subplots()
    plt.loglog(sizes, cu, 'ro-')
    plt.loglog(sizes, cv, 'bo-')
    plt.loglog([min(sizes), max(sizes)], [1e-7*min(sizes)**2, 1e-7*max(sizes)**2], 'k--')
    plt.legend(["O2", "CO2", "Quadratic convergence"])
    plt.xlabel("Element size (edge length in mm)")
    plt.ylabel("Relative L2 norm of errors")
    plt.grid()
    plt.show()

    