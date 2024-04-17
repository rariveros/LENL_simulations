import math

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.collections import PolyCollection

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    directory = "C:/mnustes_science/simulation_data/FD/localized_chaos/EE"
    PDF_EE = np.loadtxt(directory + "/PDF_EEs.txt", delimiter=',')
    dH_EE = np.loadtxt(directory + "/dH_EEs.txt", delimiter=',')
    gammas = np.loadtxt(directory + "/gammas.txt", delimiter=',')
    p_EE = np.loadtxt(directory + "/p_EE.txt", delimiter=',')

def polygon_under_graph(x, y):
    return [(0, -12), *zip(x, y), (x[-1], y[-1])]


ax = plt.figure().add_subplot(projection='3d')
mean = [np.sum(dH_EE[i] * np.exp(PDF_EE[i])) / np.sum(np.exp(PDF_EE[i])) for i in range(len(dH_EE))]
# verts[i] is a list of (x, y) pairs defining polygon i.
verts = [polygon_under_graph(dH_EE[i]/mean[i], PDF_EE[i])
         for i in range(len(dH_EE))]
facecolors = [(1 - (n / len(dH_EE)), 0, n / len(dH_EE)) for n in range(len(dH_EE))]#plt.colormaps['viridis_r'](np.linspace(0, 1, len(verts)))

poly = PolyCollection(verts, facecolors=facecolors, alpha=.7)
ax.add_collection3d(poly, zs=gammas[1:], zdir='y')

ax.set(xlim=(0, 5), ylim=(gammas[0], gammas[-1] + (gammas[1] - gammas[0])), zlim=(-12, 0),
       xlabel='$h/\langle h \\rangle$', ylabel=r'$\gamma_0$', zlabel='$\\textrm{ln[PDF]}$')
plt.savefig('EE_distibution.png', dpi=300)
plt.show()
plt.close()