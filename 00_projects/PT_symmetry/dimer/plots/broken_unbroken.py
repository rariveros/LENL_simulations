import matplotlib.pyplot as plt

from back_process import *
from functions import *
from back_process import *
from time_integrators import *

if __name__ == "__main__":
    principal_directory = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras"

    disc = "C:"
    nu = "0.018"
    gamma = "0.185"
    sigma = "6.000"
    alpha = "6.524"
    beta = 0.004811649356064012
    distances = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/dists.txt', delimiter=',')
    T = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/t_grid.txt', delimiter=',')
    x_grid = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu + '/sigma=' + sigma + '/gamma=' + gamma + '/analysis/x_grid.txt', delimiter=',')
    UR_Rs = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/UR_Rs.txt', delimiter=',')
    UR_Is = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/UR_Is.txt', delimiter=',')
    UL_Rs = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/UL_Rs.txt', delimiter=',')
    UL_Is = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/UL_Is.txt', delimiter=',')
    phi_R = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/ansatz_right.txt', delimiter=',', dtype=complex)
    phi_L = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/ansatz_left.txt', delimiter=',', dtype=complex)
    save_directory = disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis'

    for i in range(len(distances)):
        if distances[i] == 31:
            UR_R01 = UR_Rs[i]
            UR_I01 = UR_Is[i]
            UL_R01 = UL_Rs[i]
            UL_I01 = UL_Is[i]
            Z_01 = np.outer(UR_Rs[i] + 1j * UR_Is[i], phi_R[i]) + np.outer(UL_Rs[i] + 1j * UL_Is[i], phi_L[i])
        elif distances[i] == 21:
            UR_R02 = UR_Rs[i]
            UR_I02 = UR_Is[i]
            UL_R02 = UL_Rs[i]
            UL_I02 = UL_Is[i]
            Z_02 = np.outer(UR_Rs[i] + 1j * UR_Is[i], phi_R[i]) + np.outer(UL_Rs[i] + 1j * UL_Is[i], phi_L[i])
    t0 = 1000
    ti = 0
    tf = 400

    fig, ((ax3, ax1), (ax4, ax2)) = plt.subplots(nrows=2, ncols=2, figsize=(9, 4.5))

    ax1.plot(T - t0, UR_R01, color="b")
    ax1.plot(T - t0, UR_I01, color="r")
    ax1.plot(T - t0, UL_R01, color="b", linestyle="--")
    ax1.plot(T - t0, UL_I01, color="r", linestyle="--")
    ax1.set_xlim(ti, tf)
    ax1.tick_params(axis="x", direction="in", labeltop=False, labelbottom=True, top=True, bottom=True)
    ax1.tick_params(axis="y", direction="in", left=True, right=True)
    ax1.grid(alpha=0.2, color="k")

    ax2.plot(T - t0, UR_R02, color="b", label="$\\textrm{Re}\{\phi_R\}$")
    ax2.plot(T - t0, UR_I02, color="r", label="$\\textrm{Im}\{\phi_R\}$")
    ax2.plot(T - t0, UL_R02, color="b", linestyle="--", label="$\\textrm{Re}\{\phi_L\}$")
    ax2.plot(T - t0, UL_I02, color="r", linestyle="--", label="$\\textrm{Im}\{\phi_L\}$")
    ax2.set_xlim(ti, tf)
    ax2.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax2.tick_params(axis="y", direction="in", left=True, right=True)
    ax2.grid(alpha=0.2, color="k")
    ax2.legend(loc='upper left', fontsize=12)

    ax3.plot(T - t0, np.abs(UR_R01 + 1j * UR_I01)/ np.sqrt(beta), color="k")
    ax3.plot(T - t0, np.abs(UL_R01 + 1j * UL_I01)/ np.sqrt(beta), color="k", linestyle="--")
    ax3.grid(alpha=0.2, color="k")
    ax3.set_xlim(ti, tf)
    ax3.set_xlabel("$t/T$", fontsize=20)
    ax3.tick_params(axis="x", direction="in", top=True, bottom=True)
    ax3.tick_params(axis="y", direction="in", left=True, right=True)

    ax4.plot(T - t0, np.abs(UR_R02 + 1j * UR_I02)/ np.sqrt(beta), color="k", label="$|\phi_R|$")
    ax4.plot(T - t0, np.abs(UL_R02 + 1j * UL_I02)/ np.sqrt(beta), linestyle="--", color="k", label="$|\phi_L|$")
    ax4.set_xlim(ti, tf)
    ax4.grid(alpha=0.2, color="k")
    ax4.set_xlabel("$t/T$", fontsize=20)
    ax4.tick_params(axis="x", direction="in", top=True, bottom=True)
    ax4.tick_params(axis="y", direction="in", left=True, right=True, zorder=10)
    ax4.legend(loc='upper left', fontsize=12)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.2, hspace=0.07)
    plt.savefig(save_directory + "/timeeseries.png", dpi=300)
    plt.savefig("test_timeseries.png", dpi=300)
    plt.close()

    #################################
    xi_sim, xf_sim = -75, 75
    ti_sim, tf_sim = 0, 400
    xi, xf = -75, 75
    fig, ((ax4, ax2), (ax3, ax1)) = plt.subplots(nrows=2, ncols=2)

    pc_01 = ax1.pcolor(x_grid, T - t0, np.abs(Z_01)/ np.sqrt(beta), cmap=parula_map)
    ax1.set_ylim(ti, tf)
    ax1.set_xlim(xi, xf)
    ax1.set_yticklabels([])
    ax1.set_xlabel("$x\ \\textrm{(mm)}$", fontsize=20)
    ax1.grid(alpha=0.2, color="k")
    ax1.tick_params(axis="x", direction="in", labeltop=False, labelbottom=True, top=True, bottom=True)
    ax1.tick_params(axis="y", direction="in", left=True, right=True)

    ax1.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax1.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax1.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax1.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax1.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{2MA}$", fontsize=13, rotation=45, zorder=13)

    pc_02 = ax2.pcolor(x_grid, T - t0, np.abs(Z_02)/ np.sqrt(beta), cmap=parula_map)
    ax2.set_ylim(ti, tf)
    ax2.set_xlim(xi, xf)
    ax2.set_yticklabels([])
    #ax2.set_xticklabels([])
    #ax2.set_ylabel("$t/T$", fontsize=20)
    #ax2.set_xlabel("$x\ \\textrm{(mm)}$", fontsize=20)
    ax2.grid(alpha=0.2, color="k")
    ax2.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax2.tick_params(axis="y", direction="in", left=True, right=True, zorder=10)

    ax2.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax2.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k",zorder=12)
    ax2.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax2.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax2.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{2MA}$", fontsize=13, rotation=45, zorder=13)

    ax3.set_xticklabels([])
    ax3.set_yticklabels([])
    ax4.set_xticklabels([])
    ax4.set_yticklabels([])

    plt.subplots_adjust(wspace=0.2, hspace=0.07)
    cax_01 = fig.add_axes([0.913, 0.11, 0.015, 0.37])
    fig.colorbar(pc_01, cax=cax_01)
    cax_02 = fig.add_axes([0.913, 0.51, 0.015, 0.37])
    fig.colorbar(pc_02, cax=cax_02)
    plt.savefig(save_directory + "/ST_recovered.png", dpi=300)
    plt.close()