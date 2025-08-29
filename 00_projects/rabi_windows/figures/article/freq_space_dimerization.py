from functions import *
from back_process import *
from time_integrators import *

import numpy as np
import matplotlib.pyplot as plt
import os, tkinter as tk
from tkinter import filedialog

# --- parche: imports para contorno sub-pixel ---
from skimage import measure
from scipy.signal import savgol_filter  # opcional para suavizar la línea

if __name__ == '__main__':
    disc = "X:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]

    FREQS = []
    K = []
    nus = []

    for directory_01 in directories:
        print(directory_01)
        dir_01 = working_directory + "/" + directory_01 + "/sigma=3.000/gamma=0.280"
        freqs = np.loadtxt(dir_01 + '/analysis/freqs.txt', delimiter=',')
        powers = np.loadtxt(dir_01 + '/analysis/powers.txt', delimiter=',')
        nu = float(directory_01.split("=")[-1])
        FREQS.append(freqs[:, 0])
        nus.append(nu)

    dist = np.loadtxt(dir_01 + '/analysis/dists.txt', delimiter=',')
    NU, K = np.meshgrid(nus, dist)
    FREQS_array = np.transpose(np.array(FREQS))

    fig, axes = plt.subplots(1, 1, figsize=(5.5, 4))

    labelsize = 13
    xy_labelsize = 17

    pc_01a = axes.pcolor(K, NU, 1000 * FREQS_array, cmap="inferno", vmin=0, vmax= 1000 * np.amax(FREQS_array[:, :-9]))
    cax_01a = fig.add_axes([0.87, 0.2, 0.025, 0.6])
    cb1a = fig.colorbar(pc_01a, cax=cax_01a)
    cb1a.set_label('$\Omega \\times 10^{-3}$', rotation=0, size=xy_labelsize, labelpad=-30, y=1.2)
    cb1a.ax.tick_params(labelsize=labelsize)

    level = 7*1e-4
    contours = measure.find_contours(FREQS_array, level=level)  # devuelve curvas en coords (row, col)
    dist_arr = np.asarray(dist)
    nus_arr = np.asarray(nus)

    for C in contours:
        rows, cols = C[:, 0], C[:, 1]
        K_line  = np.interp(rows, np.arange(len(dist_arr)), dist_arr)
        NU_line = np.interp(cols, np.arange(len(nus_arr)),  nus_arr)

        if len(K_line) >= 31:
            K_line  = savgol_filter(K_line, 15, 3)
            NU_line = savgol_filter(NU_line, 20, 3)

        axes.plot(K_line, NU_line, lw=4.5, color='w', zorder=10)
        axes.plot(K_line, NU_line, lw=1.5, color='k', zorder=11)
    # ---------------------------------------------------------------------
    axes.tick_params(axis="x", direction="in", labelsize=labelsize, top=True, bottom=True, labeltop=True, labelbottom=False)
    axes.tick_params(axis="y", direction="in", labelsize=labelsize, left=True, right=True, labelleft=True, labelright=False)
    axes.set_xticks([10, 20, 30, 40])
    axes.set_ylim(0, 0.32)
    axes.set_xlabel("$d$", fontsize=xy_labelsize)
    axes.set_ylabel("$\\nu$", fontsize=xy_labelsize)
    axes.grid(alpha=0.4, ls="--")
    fig.subplots_adjust(wspace=0.35, hspace=0.25, left=0.15, right=0.85, bottom=0.2, top=0.8)
    plt.savefig("parameters_space_dimerization.png", dpi=300)
    plt.close()