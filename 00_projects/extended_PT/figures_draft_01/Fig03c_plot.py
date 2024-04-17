from back_process import *

if __name__ == "__main__":
    directory = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis"
    C_r = np.loadtxt(directory + '/coef_real_dist.txt', delimiter=',')
    C_i = np.loadtxt(directory + '/coef_imag_dist.txt', delimiter=',')
    dist = np.loadtxt(directory + '/coef_dist.txt', delimiter=',')

    fig, ((ax2), (ax1), (ax3)) = plt.subplots(ncols=1, nrows=3)
    ax1.set_xlabel('$d\ \\textrm{(mm)}$', fontsize=25)
    ax1.set_ylabel('$\\xi$', fontsize=25)
    ax1.tick_params(labelsize=15)
    ax1.scatter(dist, C_r, zorder=10, label="$\\textrm{Re}(\\xi)$")
    ax1.scatter(dist, C_i, zorder=10, label="$\\textrm{Im}(\\xi)$")
    ax1.grid(alpha=0.2, zorder=0)
    ax1.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig("Fig03c.png", dpi=300)
    plt.close()