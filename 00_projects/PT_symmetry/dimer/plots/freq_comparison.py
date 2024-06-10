from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    save_directory_01a = 'C:/mnustes_science/simulation_data/FD/PT_dimer/alpha=6.524/beta=1.000/mu=0.100/nu=0.020/sigma=6.000/gamma=0.185/analysis'
    directory = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis"

    #data = np.loadtxt(directory + '/data.txt', delimiter=',')
    data_dist = np.loadtxt(directory + '/data_dist.txt', delimiter=',')

    init = 0
    fin = -1

    ### FREQUENCY VS DISTANCE ###
    dist_full = data_dist[:fin, 0]
    freq_dist_full = 2 * np.pi * 1000 / data_dist[:fin, 1]
    print(dist_full)
    print(freq_dist_full)
    freq_dist_err_full = 2 * np.pi * 1000 * np.abs(data_dist[:fin, 2]/data_dist[:fin, 1]**2)

    init = -20
    dist_bif = dist_full[init:]
    freq_dist_bif = freq_dist_full[init:]
    freq_dist_err_bif = freq_dist_err_full[init:]

    dist_01 = dist_full[1:21]
    freq_dist_01 = freq_dist_full[1:21]
    freq_dist_err_01 = freq_dist_err_full[1:21]

    dist_02 = dist_full[0:46]
    freq_dist_02 = freq_dist_full[0:46]
    freq_dist_err_02 = freq_dist_err_full[0:46]


    delta_d = dist_bif[1] - dist_bif[0]

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('$d\ \\textrm{(mm)}$', fontsize=25)
    ax1.set_ylabel('$\Omega \\times 10^{-3}$', fontsize=25)
    ax1.tick_params(labelsize=20)
    ax1.scatter(dist_full, freq_dist_full, c="k", s=40, zorder=10, label="$\\textbf{NUM}$")
    zeros_01 = np.arange(30.3, 35, 0.25)
    zeros_02 = np.arange(14, 18.5, 0.25)
    ax1.scatter(zeros_01, 0 * zeros_01, c="k", s=20, zorder=10)
    ax1.scatter(zeros_02, 0 * zeros_02, c="k", s=20, zorder=10)
    ax1.hlines(0, 0, 33, colors="k")
    ax1.vlines(18.41, 0, 150, colors="k", linestyles="--")
    ax1.vlines(30.2930, 0, 150, colors="k", linestyles="--")
    ax1.set_xlim([16, 32.5])
    ax1.set_ylim([-1, 110])
    ax1.grid(alpha=0.2, zorder=0)

    distances_01a = np.loadtxt(save_directory_01a + '/dists.txt', delimiter=',')
    freqs_01a = np.loadtxt(save_directory_01a + '/freqs.txt', delimiter=',')
    alpha = 1

    ax1.plot(distances_01a[2:], 2 * freqs_01a[2:, 0], linewidth=2, color="r", linestyle="--", label="$\\textbf{2MA}$")
    ax1.legend(fontsize=20)
    plt.tight_layout()
    #plt.savefig("freqs.png", dpi=150)
    plt.savefig("freq_comparisson.png", dpi=200)
    plt.close()