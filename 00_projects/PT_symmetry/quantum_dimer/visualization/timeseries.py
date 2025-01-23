from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    ############ DAMPED ############
    working_directory = "D:/mnustes_science/simulation_data/FD/ladder_operators/test/Delta=0.1000/gamma=0.1000"

    U1_light = np.loadtxt(working_directory + '/Omega=0.0500/k=0.1400/U1.txt', delimiter=',', dtype=np.complex128)
    V1_light = np.loadtxt(working_directory + '/Omega=0.0500/k=0.1400/V1.txt', delimiter=',', dtype=np.complex128)
    t_light = np.loadtxt(working_directory + '/Omega=0.0500/k=0.1400/T.txt', delimiter=',')
    ti = 0
    tf = 500
    fig, (ax1, ax2) = plt.subplots(2, figsize=(6, 2.5))
    ax1.plot(t_light[:-1], np.real(U1_light) * 1000, c="#5489ef", lw=3, zorder=4)
    ax1.plot(t_light[:-1], np.imag(U1_light) * 1000, c="#0900bc", lw=3, zorder=5)
    ax1.hlines(0, 0, 500, colors="k", alpha=1.0, zorder=2)
    ax1.tick_params("both", direction="in", labelsize=18, labelbottom=False, labeltop=False, labelleft=True, labelright=False)
    ax1.set_ylabel('$\\beta_1$', size=25)
    ax1.set_xlim([ti, tf])
    ax1.set_ylim([-10, 10])
    ax1.set_yticks([-8.0, 0.0, 8.0], ["$-8.0$", "$0.0$", "$8.0$"])
    ax1.grid(linestyle='--', alpha=0.5)
    ax1

    ax2.plot(t_light[:-1], np.real(V1_light) * 1000, c="#ff6767", lw=3, zorder=5)
    ax2.plot(t_light[:-1], np.imag(V1_light) * 1000, c="#980012", lw=3, zorder=4)
    ax2.hlines(0, 0, 500, colors="k", alpha=1.0, zorder=2)
    ax2.tick_params("both", direction="in", labelsize=18, labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    ax2.set_xlabel('$t$', size=25)
    ax2.set_ylabel('$\\beta_2$', size=25)
    ax2.set_xlim([ti, tf])
    ax2.set_ylim([-10, 10])
    ax2.set_yticks([-8.0, 0.0, 8.0], ["$-8.0$", "$0.0$", "$8.0$"])
    ax2.grid(linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig("timeseries_damped.png", dpi=300)
    plt.close()

    ############ RABI ############

    working_directory = "D:/mnustes_science/simulation_data/FD/ladder_operators/test/Delta=0.1000/gamma=0.1000"

    U1_light = np.loadtxt(working_directory + '/Omega=0.0800/k=0.1400/U1.txt', delimiter=',', dtype=np.complex128)
    V1_light = np.loadtxt(working_directory + '/Omega=0.0800/k=0.1400/V1.txt', delimiter=',', dtype=np.complex128)
    t_light = np.loadtxt(working_directory + '/Omega=0.0800/k=0.1400/T.txt', delimiter=',')
    ti = 0
    tf = 500
    fig, (ax1, ax2) = plt.subplots(2, figsize=(6, 2.5))
    ax1.plot(t_light[:-1], np.real(U1_light), c="#5489ef", lw=3, zorder=4)
    ax1.plot(t_light[:-1], np.imag(U1_light), c="#0900bc", lw=3, zorder=5)
    ax1.hlines(0, 0, 500, colors="k", alpha=1.0, zorder=2)
    ax1.tick_params("both", direction="in", labelsize=18, labelbottom=False, labeltop=False, labelleft=True, labelright=False)
    ax1.set_ylabel('$\\beta_1$', size=25)
    ax1.set_xlim([ti, tf])
    ax1.set_ylim([-0.25, 0.25])
    ax1.set_yticks([-0.2, 0, 0.2])
    ax1.grid(linestyle='--', alpha=0.5)
    ax1

    ax2.plot(t_light[:-1], np.real(V1_light), c="#ff6767", lw=3, zorder=5)
    ax2.plot(t_light[:-1], np.imag(V1_light), c="#980012", lw=3, zorder=4)
    ax2.hlines(0, 0, 500, colors="k", alpha=1.0, zorder=2)
    ax2.tick_params("both", direction="in", labelsize=18, labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    ax2.set_xlabel('$t$', size=25)
    ax2.set_ylabel('$\\beta_2$', size=25)
    ax2.set_xlim([ti, tf])
    ax2.set_ylim([-0.25, 0.25])
    ax2.set_yticks([-0.2, 0, 0.2])
    ax2.grid(linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig("timeseries_rabi.png", dpi=300)
    plt.close()

    ############ STATIONARY ############

    working_directory = "D:/mnustes_science/simulation_data/FD/ladder_operators/test/Delta=0.1000/gamma=0.1000"

    U1_light = np.loadtxt(working_directory + '/Omega=0.0800/k=0.0000/U1.txt', delimiter=',', dtype=np.complex128)
    V1_light = np.loadtxt(working_directory + '/Omega=0.0800/k=0.0000/V1.txt', delimiter=',', dtype=np.complex128)
    t_light = np.loadtxt(working_directory + '/Omega=0.0800/k=0.0000/T.txt', delimiter=',')
    ti = 0
    tf = 500
    fig, (ax1, ax2) = plt.subplots(2, figsize=(6, 2.5))
    ax1.plot(t_light[:-1], np.real(U1_light), c="#5489ef", lw=3, zorder=4)
    ax1.plot(t_light[:-1], np.imag(U1_light), c="#0900bc", lw=3, zorder=5)
    ax1.hlines(0, 0, 500, colors="k", alpha=1.0, zorder=2)
    ax1.tick_params("both", direction="in", labelsize=18, labelbottom=False, labeltop=False, labelleft=True, labelright=False)
    ax1.set_ylim([-0.25, 0.25])
    ax1.set_yticks([-0.2, 0, 0.2])
    ax1.set_ylabel('$\\beta_1$', size=25)
    ax1.set_xlim([ti, tf])
    ax1.grid(linestyle='--', alpha=0.5)

    ax2.plot(t_light[:-1], np.real(V1_light), c="#ff6767", lw=3, zorder=5)
    ax2.plot(t_light[:-1], np.imag(V1_light), c="#980012", lw=3, zorder=4)
    ax2.hlines(0, 0, 500, colors="k", alpha=1.0, zorder=2)
    ax2.tick_params("both", direction="in", labelsize=18, labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    ax2.set_yticks([-0.2, 0 , 0.2])
    ax2.set_xlabel('$t$', size=25)
    ax2.set_ylabel('$\\beta_2$', size=25)
    ax2.set_xlim([ti, tf])
    ax2.set_ylim([-0.25, 0.25])
    ax2.grid(linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig("timeseries_stationary.png", dpi=300)
    plt.close()