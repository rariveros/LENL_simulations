import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':
    f1, f2 = 0.480, 0.530
    T = 1 / np.abs((f1 - f2))
    t_grid = np.arange(0, T, 0.01)
    x1 = np.cos(2 * np.pi * f1 * (t_grid - 0.25 * T))
    x2 = np.cos(2 * np.pi * f2 * (t_grid - 0.25 * T))
    x3 = x1 + x2

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7, 5))
    ax1.plot(t_grid, x1, color='k', lw=2)
    ax2.plot(t_grid, x2, color='k', lw=2)
    ax3.plot(t_grid, x3, color='k', lw=2)
    ax1.grid(axis='y')
    ax2.grid(axis='y')
    ax3.grid(axis='y')
    ax1.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=False, labelright=False)
    ax1.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax2.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=False, labelright=False)
    ax2.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax3.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=False, labelright=False)
    ax3.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax1.set_xlim(t_grid[0], t_grid[-1])
    ax2.set_xlim(t_grid[0], t_grid[-1])
    ax3.set_xlim(t_grid[0], t_grid[-1])
    plt.savefig("wave_profile.png", dpi=300)