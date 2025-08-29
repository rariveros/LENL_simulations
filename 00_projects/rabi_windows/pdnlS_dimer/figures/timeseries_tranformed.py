from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    dir_RO = r"C:\mnustes_science\simulation_data\FD\pdnlS_dimer\nu=0.1000\mu=0.1000\gamma=0.2000\k=0.0200"
    # dir_RO = r"C:\mnustes_science\simulation_data\FD\pdnlS_dimer\nu=0.1000\mu=0.1000\gamma=0.2000\k=0.0700"
    # dir_steady = r"C:\mnustes_science\simulation_data\FD\pdnlS_dimer\nu=0.1000\mu=0.1000\gamma=0.2000\k=0.0200"
    # dir_damped = r"C:\mnustes_science\simulation_data\FD\pdnlS_dimer\nu=0.2000\mu=0.1000\gamma=0.2000\k=0.0200"

    T_RO = np.loadtxt(dir_RO + '/T.txt', delimiter=',')[:-1]
    Z1_RO = np.loadtxt(dir_RO + '/U.txt', delimiter=',', dtype=np.complex128)
    Z2_RO = np.loadtxt(dir_RO + '/V.txt', delimiter=',', dtype=np.complex128)

    U_RO = 0.5 * (Z1_RO + Z2_RO + np.conjugate(Z1_RO - Z2_RO))
    V_RO = 0.5 * (Z1_RO + Z2_RO - np.conjugate(Z1_RO - Z2_RO))

    #T_steady = np.loadtxt(dir_steady + '/T.txt', delimiter=',')[:-1]
    #U_steady = np.loadtxt(dir_steady + '/U.txt', delimiter=',', dtype=np.complex128)
    #V_steady = np.loadtxt(dir_steady + '/V.txt', delimiter=',', dtype=np.complex128)


    t0 = 100
    fig, (ax1, ax2) = plt.subplots(2,figsize=(3., 2))
    ax1.plot(T_RO - t0, np.real(U_RO), c="#6e04bf", label="$\\textrm{Re }z_1$", lw=2, zorder=5)
    ax1.plot(T_RO - t0, np.imag(U_RO), c="#55d604", label="$\\textrm{Im }z_1$", lw=2, zorder=5)
    ax1.plot(T_RO - t0, np.abs(U_RO), c="k", label="$|z_1|$", lw=2, zorder=5)
    ax1.hlines(0, 0, 400, colors="k", linewidths=0.7)
    ax1.set_ylabel('$u$', size=16)
    ax1.set_xlim(0, 400)
    ax1.tick_params(axis="both", direction="in", labelsize=13, labelbottom=False)
    ax1.set_yticks([-0.2, 0, 0.2])
    ax1.set_ylim(-0.35, 0.35)
    #ax1.legend(fontsize=16)
    ax1.grid(alpha=0.3)

    ax2.plot(T_RO - t0, np.real(V_RO), c="#6e04bf", label="$\\textrm{Re }z_2$", lw=2, zorder=5)
    ax2.plot(T_RO - t0, np.imag(V_RO), c="#55d604", label="$\\textrm{Im }z_2$", lw=2, zorder=5)
    ax2.plot(T_RO - t0, np.abs(V_RO), c="k", label="$|z_2|$", lw=2, zorder=5)
    ax2.hlines(0, 0, 400, colors="k", linewidths=0.7)
    ax2.set_ylabel('$v$', size=16)
    ax2.set_xlabel('$\\textrm{Time}$', size=16)
    ax2.set_xlim(0, 400)
    ax2.set_ylim(-0.35, 0.35)
    ax2.tick_params(axis="both", direction="in", labelsize=13)
    ax2.set_yticks([-0.2, 0, 0.2])
    ax2.grid(alpha=0.3)

    fig.subplots_adjust(left=0.3, right=0.95, bottom=0.3, top=0.8)
    plt.savefig("timeseries_transformed.png", dpi=300)
    plt.close()