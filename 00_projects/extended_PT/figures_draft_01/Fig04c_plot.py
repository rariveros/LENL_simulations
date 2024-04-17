from back_process import *

if __name__ == "__main__":
    directory = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis"
    C_r = np.loadtxt(directory + '/coef_real_gamma.txt', delimiter=',')
    C_i = np.loadtxt(directory + '/coef_imag_gamma.txt', delimiter=',')
    dist = np.loadtxt(directory + '/coef_gamma.txt', delimiter=',')

    fig, (ax1) = plt.subplots(ncols=1, nrows=1)
    ax1.set_xlabel('$\gamma_0$', fontsize=35)
    ax1.set_ylabel('$\langle \\xi \\rangle$', fontsize=35)
    ax1.tick_params(labelsize=25)
    ax1.scatter(dist, C_r, zorder=10, label="$\\textrm{Re}(\\xi)$")
    ax1.scatter(dist, C_i, zorder=10, label="$\\textrm{Im}(\\xi)$")
    ax1.grid(alpha=0.2, zorder=0)
    ax1.legend(fontsize=20)
    ax1.set_aspect(0.06)
    plt.tight_layout()
    plt.savefig("Fig04c.png", dpi=300)
    plt.close()