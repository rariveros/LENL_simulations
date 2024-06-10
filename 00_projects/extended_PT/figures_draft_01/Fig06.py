
import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    directory = "C:/Users/research/LENL_simulations/00_projects/PT_symmetry/spectral_analysis/spectral_data"
    beta = 0.004811649356064012

    mode_evol_zero = np.loadtxt(directory + '/mod_ev_zero.txt', delimiter=',')
    mode_0_zero = np.loadtxt(directory + '/mod_0_zero.txt', delimiter=',')
    eig_r_zero = np.loadtxt(directory + '/eig_r_zero.txt', delimiter=',')
    eig_i_zero = np.loadtxt(directory + '/eig_i_zero.txt', delimiter=',')
    DIST_zero = np.loadtxt(directory + '/DIST_zero.txt', delimiter=',')
    x_grid_zero = np.loadtxt(directory + '/x_grid_zero.txt', delimiter=',')

    mode_evol_sol = np.loadtxt(directory + '/mod_ev_sol.txt', delimiter=',')
    mode_0_sol = np.loadtxt(directory + '/mod_0_sol.txt', delimiter=',')
    eig_r_sol = np.loadtxt(directory + '/eig_r_sol.txt', delimiter=',')
    eig_i_sol = np.loadtxt(directory + '/eig_i_sol.txt', delimiter=',')
    DIST_sol = np.loadtxt(directory + '/DIST_sol.txt', delimiter=',')
    x_grid_sol = np.loadtxt(directory + '/x_grid_sol.txt', delimiter=',')

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2,gridspec_kw={'width_ratios': [1.5, 1.0]})
    plt.subplots_adjust(wspace=0.08, hspace=0.1)

    n = 100

    sc = ax1.scatter(eig_r_zero * n, eig_i_zero* n, edgecolors="0.", c=DIST_zero, zorder=10, cmap=parula_map)
    ax1.scatter(np.flip(eig_r_zero)[6] * n, np.flip(eig_i_zero)[6] * n, linewidth=1.5, edgecolors="0.", c="r", zorder=12, s=85)
    ax1.arrow(-0.0011 * n, 0.0365 * n, 0.002 * n, 0, width=0.003 * n, head_length=0.001 * n, head_width=0.01 * n, facecolor="r", zorder=11, edgecolor="0.")
    ax1.arrow(-0.0011 * n, -0.0365 * n, 0.002 * n, 0, width=0.003 * n, head_length=0.001 * n, head_width=0.01 * n, facecolor="r", zorder=11, edgecolor="0.")
    cbar_01 = plt.colorbar(sc)
    cbar_01.set_label('$d\ \\textrm{(mm)}$', rotation=0, size=15, labelpad=-10, y=1.15)
    cbar_01.ax.tick_params(labelsize=11)
    ax1.set_xlim([-0.015 * n, 0.003 * n])
    ax1.set_ylim([-0.06 * n, 0.06 * n])
    ax1.vlines(0, -20, 20, colors="k", alpha=0.8)
    ax1.hlines(0, -20, 20, colors="k", alpha=0.8)
    ax1.grid(alpha=0.4, zorder=0)
    ax1.set_ylabel("$\\textrm{Im}\{\lambda_i\}\\times 10^{-2}$", fontsize=15)
    ax1.tick_params(axis="x", direction="in", top=True, bottom=True, labeltop=True, labelbottom=False, labelsize=12)
    ax1.tick_params(axis="y", direction="in", left=True, right=True, labelleft=True, labelright=False, labelsize=12)

    ax2.plot(x_grid_zero, mode_evol_zero / np.sqrt(beta), color="k", label="$|A_s + e^{\lambda_c dt}A_c|$", linewidth=2)
    ax2.plot(x_grid_zero, mode_0_zero / np.sqrt(beta), color="k", label="$|A_s|$", linestyle=":", linewidth=2)
    ax2.set_xlim(-75, 75)
    ax2.grid(alpha=0.3)
    ax2.tick_params(axis="x", direction="in", top=True, bottom=True, labeltop=True, labelbottom=False, labelsize=12)
    ax2.tick_params(axis="y", direction="in", left=True, right=True, labelleft=False, labelright=True, labelsize=12)
    ax2.set_ylabel("$|A|\ \\textrm{(mm)}$", fontsize=15, rotation=270, labelpad=20)
    ax2.yaxis.set_label_position("right")

    sc = ax3.scatter(np.flip(eig_r_sol) * n, np.flip(eig_i_sol) * n, edgecolors="0.", c=np.flip(DIST_sol), zorder=10, cmap=parula_map)
    ax3.scatter(np.flip(eig_r_sol)[-1] * n, np.flip(eig_i_sol)[-1] * n, linewidth=1.5, edgecolors="0.", c="r", zorder=12, s=85)
    ax3.arrow(-0.006 * n, 0 * n, 0.009 * n, 0. * n, width=0.003 * n, head_length=0.005 * n, head_width=0.01 * n, facecolor="r", zorder=11, edgecolor="0.")
    cbar_03 = plt.colorbar(sc)
    cbar_03.ax.tick_params(labelsize=11)
    #cbar.set_label('$d\ \\textrm{(mm)}$', rotation=0, size=15, labelpad=20, y=1.1)
    ax3.set_xlim([-0.07 * n, 0.014 * n])
    ax3.set_ylim([-0.06 * n, 0.06 * n])
    ax3.vlines(0, -20, 20, colors="k", alpha=0.8)
    ax3.hlines(0, -20, 20, colors="k", alpha=0.8)
    ax3.grid(alpha=0.4, zorder=0)
    ax3.set_xlabel("$\\textrm{Re}\{\lambda_i\}\\times 10^{-2}$", fontsize=15)
    ax3.set_ylabel("$\\textrm{Im}\{\lambda_i\}\\times 10^{-2}$", fontsize=15)
    ax3.tick_params(axis="x", direction="in", top=True, bottom=True, labelsize=12)
    ax3.tick_params(axis="y", direction="in", left=True, right=True, labelsize=12)

    ax4.plot(x_grid_sol, mode_evol_sol / np.sqrt(beta), color="k", label="$|A_s + e^{\lambda_c dt}A_c|$", linewidth=2)
    ax4.plot(x_grid_sol, mode_0_sol / np.sqrt(beta), color="k", label="$|A_s|$", linestyle=":", linewidth=2)
    ax4.set_xlim(-75, 75)
    ax4.grid(alpha=0.3)
    ax4.set_xlabel("$x\ \\textrm{(mm)}$", fontsize=15)
    ax4.tick_params(axis="x", direction="in", top=True, bottom=True, labelsize=12, labeltop=False, labelbottom=True)
    ax4.tick_params(axis="y", direction="in", left=True, right=True, labelsize=12, labelleft=False, labelright=True)
    ax4.set_ylabel("$|A|\ \\textrm{(mm)}$", fontsize=15, rotation=270, labelpad=20)
    ax4.yaxis.set_label_position("right")

    props = dict(boxstyle='round,pad=0.1', facecolor='white', alpha=0.75)
    ax1.text(-1.45, 4.9, "$\\textbf{(a)}$", fontsize=14, rotation=0, zorder=13)
    ax2.text(-70.3, 4.7, "$\\textbf{(b)}$", fontsize=14, rotation=0, zorder=13)
    ax3.text(-6.77, 4.9, "$\\textbf{(c)}$", fontsize=14, rotation=0, zorder=13)
    ax4.text(-70.3, 5, "$\\textbf{(d)}$", fontsize=14, rotation=0, zorder=13)

    #savedirectory = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis/figures"
    plt.savefig('Fig06.png', dpi=300)
    plt.close()
