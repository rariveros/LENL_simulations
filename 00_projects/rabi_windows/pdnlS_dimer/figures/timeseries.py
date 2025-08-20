from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    dir_RO = r"C:\mnustes_science\simulation_data\FD\pdnlS_dimer\nu=0.1000\mu=0.1000\gamma=0.2000\k=0.0700"
    # dir_RO = r"C:\mnustes_science\simulation_data\FD\pdnlS_dimer\nu=0.1000\mu=0.1000\gamma=0.2000\k=0.0700"
    # dir_steady = r"C:\mnustes_science\simulation_data\FD\pdnlS_dimer\nu=0.1000\mu=0.1000\gamma=0.2000\k=0.0200"
    # dir_damped = r"C:\mnustes_science\simulation_data\FD\pdnlS_dimer\nu=0.2000\mu=0.1000\gamma=0.2000\k=0.0200"

    T_RO = np.loadtxt(dir_RO + '/T.txt', delimiter=',')[:-1]
    U_RO = np.loadtxt(dir_RO + '/U.txt', delimiter=',', dtype=np.complex128)
    V_RO = np.loadtxt(dir_RO + '/V.txt', delimiter=',', dtype=np.complex128)

    #T_steady = np.loadtxt(dir_steady + '/T.txt', delimiter=',')[:-1]
    #U_steady = np.loadtxt(dir_steady + '/U.txt', delimiter=',', dtype=np.complex128)
    #V_steady = np.loadtxt(dir_steady + '/V.txt', delimiter=',', dtype=np.complex128)

    fig, (ax1, ax2) = plt.subplots(2,figsize=(4, 2.5))
    ax1.plot(T_RO, np.real(U_RO), c="b", label="$\\textrm{Re }z_1$", lw=3, zorder=5)
    ax1.plot(T_RO, np.imag(U_RO), c="r", label="$\\textrm{Im }z_1$", lw=3, zorder=5)
    ax1.plot(T_RO, np.abs(U_RO), c="k", label="$|z_1|$", lw=3, zorder=5)
    ax1.hlines(0, 0, 300, colors="k")
    ax1.set_ylabel('$z_1$', size=20)
    ax1.set_xlim(0, 300)
    ax1.tick_params(axis="both", direction="in", labelsize=16, labelbottom=False)
    ax1.set_yticks([-0.2, 0, 0.2])
    #ax1.legend(fontsize=16)
    ax1.grid(alpha=0.3)

    ax2.plot(T_RO, np.real(V_RO), c="b", label="$\\textrm{Re }z_2$", lw=3, zorder=5)
    ax2.plot(T_RO, np.imag(V_RO), c="r", label="$\\textrm{Im }z_2$", lw=3, zorder=5)
    ax2.plot(T_RO, np.abs(V_RO), c="k", label="$|z_2|$", lw=3, zorder=5)
    ax2.hlines(0, 0, 300, colors="k")
    ax2.set_ylabel('$z_2$', size=20)
    ax2.set_xlabel('$t$', size=20)
    ax2.set_xlim(0, 300)
    ax2.tick_params(axis="both", direction="in", labelsize=16)
    ax2.set_yticks([-0.2, 0, 0.2])
    ax2.grid(alpha=0.3)

    fig.subplots_adjust(left=0.22, right=0.95, bottom=0.22, top=0.98)
    plt.savefig("timeseries.png", dpi=300)
    plt.close()