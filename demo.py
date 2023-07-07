#!/usr/bin/env python3
import matplotlib as mpl
from matplotlib import tri
import matplotlib.pyplot as plt
import time
import numpy as np
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import subprocess
import sys
import os


def plot_concentrations(tr, c):
    M = int(len(c)/2)
    fig, (axu,axv) = plt.subplots(1, 2, sharex=True, sharey=True)
    
    axu.set_title('O2 (mol/m³)')
    axv.set_title('CO2 (mol/m³)')
    axu.set_xlabel('r (m)')
    axu.set_ylabel('z (m)')
    axv.set_xlabel('r (m)')
    axv.set_ylabel('z (m)')

    tpcu = axu.tripcolor(tr, c[:M], shading="gouraud")
    tpcv = axv.tripcolor(tr, c[M:], shading="gouraud")
    axu.triplot(tr, color="black", linewidth="0.1")
    axv.triplot(tr, color="black", linewidth="0.1")
    fig.colorbar(tpcu, ax=axu)
    fig.colorbar(tpcv, ax=axv)
    plt.show()

if __name__ == "__main__":
    mesh_arg = "uniform_1mm"
    constants_arg = "orchard"
    args = sys.argv
    if len(args) >= 2:
        mesh_arg = args[1]
        if len(args) == 3:
            constants_arg = args[2]

    os.chdir('C++')
    p = subprocess.run(f"./demo {mesh_arg} {constants_arg}", shell=True)

    os.chdir('..')
    points_file = f"matlab/pear_meshes/data/points_{mesh_arg}.txt"
    triangles_file = f"matlab/pear_meshes/data/triangles_{mesh_arg}.txt"
    
    rz = np.loadtxt(points_file, delimiter=',')
    r = rz[:, 0]
    z = rz[:, 1]
    triangles = np.loadtxt(triangles_file, delimiter=',', dtype=int)
    triangles -= 1
    c = np.loadtxt("c_solution.txt", delimiter=',')
    
    M = int(len(c)/2)
    tr = tri.Triangulation(r, z, triangles)
    plot_concentrations(tr, c)
    
    

