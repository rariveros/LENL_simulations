from functions import *
from back_process import *
from time_integrators import *

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
import matplotlib as mpl

if __name__ == '__main__':
    # === cargar datos ===
    dirpath = r"C:\mnustes_science\simulation_data\FD\coherent\parameter_space_01"
    output_data = np.loadtxt(dirpath + r"\output_data.txt", delimiter=',')

    # columnas: [k, quantumness, g, f_peak, M]
    k_vals = output_data[:, 0]
    q_vals = 0.1 / output_data[:, 2]
    f_peak = np.abs(output_data[:, 3])
    S_peak = np.abs(output_data[:, 4])

    # asumo que los datos están ordenados como loops (q externo, k interno)
    k_unique = np.unique(k_vals)
    q_unique = np.flip(np.unique(q_vals))

    # reshape para alinear ejes (sin doble transpose)
    f_grid = f_peak.reshape(len(q_unique), len(k_unique))
    S_grid = S_peak.reshape(len(q_unique), len(k_unique))

    S_grid = filtro_superficie(S_grid, 2, "XY")
    f_grid = filtro_superficie(f_grid, 2, "XY")

    # plot
    fig, (ax01, ax02) = plt.subplots(2, 1, figsize=(4.0, 6.5))
    pc01 = ax01.pcolormesh(k_unique, q_unique, f_grid, shading='auto', cmap='inferno', vmin= 0.0, vmax=0.18)
    ax01.tick_params(axis="both", direction="in", labelsize=11, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax01.set_ylabel(r"$\log_{10}(\gamma/g)$", fontsize=15)
    ax01.set_yscale('log')
    pc02 = ax02.pcolormesh(k_unique, q_unique, S_grid, shading='auto', cmap='inferno', norm=LogNorm(vmin=1e8, vmax=1e12))
    ax02.set_xlabel(r"$\kappa$", fontsize=15)
    ax02.set_ylabel(r"$\log_{10}(\gamma/g)$", fontsize=15)

    cax01 = fig.add_axes([0.82, 0.526, 0.025, 0.354])
    cb01 = fig.colorbar(pc01, cax=cax01)
    cb01.set_label(r"$\omega$", fontsize=14)
    cb01.ax.tick_params(labelsize=10)

    cax02 = fig.add_axes([0.82, 0.11, 0.025, 0.354])
    cb02 = fig.colorbar(pc02, cax=cax02)
    cb02.set_label(r"$S_{\textrm{max}}$", fontsize=14)
    ax02.set_yscale('log')
    cb02.ax.tick_params(labelsize=10)

    fig.subplots_adjust(wspace=0.15, hspace=0.18, left=0.26, right=0.8)
    #plt.tight_layout()
    plt.savefig(r"parameter_space_map.png", dpi=300)
    #plt.show()