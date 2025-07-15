from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    save_directory_01a = 'C:/mnustes_science/simulation_data/FD/PT_dimer/alpha=6.524/beta=1.000/mu=0.100/nu=0.020/sigma=6.000/gamma=0.185/analysis'
    save_directory_01b = 'C:/mnustes_science/simulation_data/FD/PT_dimer/alpha=6.524/beta=1.000/mu=0.100/nu=0.020/sigma=9.000/gamma=0.185/analysis'
    save_directory_01c = 'C:/mnustes_science/simulation_data/FD/PT_dimer/alpha=6.524/beta=1.000/mu=0.100/nu=0.020/sigma=12.000/gamma=0.185/analysis'
    save_directory_02 = 'C:/mnustes_science/simulation_data/FD/PT_dimer/alpha=6.524/beta=1.000/mu=0.100/nu=0.025/sigma=6.000/gamma=0.185/analysis'
    save_directory_03 = 'C:/mnustes_science/simulation_data/FD/PT_dimer/alpha=6.524/beta=1.000/mu=0.100/nu=0.030/sigma=6.000/gamma=0.185/analysis'
    save_directory_04 = 'C:/mnustes_science/simulation_data/FD/PT_dimer/alpha=6.524/beta=1.000/mu=0.100/nu=0.035/sigma=6.000/gamma=0.185/analysis'
    save_directory_05 = 'C:/mnustes_science/simulation_data/FD/PT_dimer/alpha=6.524/beta=1.000/mu=0.100/nu=0.040/sigma=6.000/gamma=0.185/analysis'
    save_directory_06 = 'C:/mnustes_science/simulation_data/FD/PT_dimer/alpha=6.524/beta=1.000/mu=0.100/nu=0.045/sigma=6.000/gamma=0.185/analysis'
    save_directory_07 = 'C:/mnustes_science/simulation_data/FD/PT_dimer/alpha=6.524/beta=1.000/mu=0.100/nu=0.050/sigma=6.000/gamma=0.185/analysis'

    save_directory_1000 = 'C:/mnustes_science/simulation_data/FD/PT_dimer/alpha=6.524/beta=1.000/mu=0.100/nu=0.100/sigma=16.000/gamma=0.240/analysis'

    distances_01a = np.loadtxt(save_directory_01a + '/dists.txt', delimiter=',')
    freqs_01a = np.loadtxt(save_directory_01a + '/freqs.txt', delimiter=',')
    #distances_02 = np.loadtxt(save_directory_02 + '/dists.txt', delimiter=',')
    #freqs_02 = np.loadtxt(save_directory_02 + '/freqs.txt', delimiter=',')
    distances_03 = np.loadtxt(save_directory_03 + '/dists.txt', delimiter=',')
    freqs_03 = np.loadtxt(save_directory_03 + '/freqs.txt', delimiter=',')
    #distances_04 = np.loadtxt(save_directory_04 + '/dists.txt', delimiter=',')
    #freqs_04 = np.loadtxt(save_directory_04 + '/freqs.txt', delimiter=',')
    distances_05 = np.loadtxt(save_directory_05 + '/dists.txt', delimiter=',')
    freqs_05 = np.loadtxt(save_directory_05 + '/freqs.txt', delimiter=',')
    #distances_06 = np.loadtxt(save_directory_06 + '/dists.txt', delimiter=',')
    #freqs_06 = np.loadtxt(save_directory_06 + '/freqs.txt', delimiter=',')
    distances_07 = np.loadtxt(save_directory_07 + '/dists.txt', delimiter=',')
    freqs_07 = np.loadtxt(save_directory_07 + '/freqs.txt', delimiter=',')

    #distances_01b = np.loadtxt(save_directory_01b + '/dists.txt', delimiter=',')
    #freqs_01b = np.loadtxt(save_directory_01b + '/freqs.txt', delimiter=',')
    #distances_01c = np.loadtxt(save_directory_01c + '/dists.txt', delimiter=',')
    #freqs_01c = np.loadtxt(save_directory_01c + '/freqs.txt', delimiter=',')
    #distances_1000 = np.loadtxt(save_directory_1000 + '/dists.txt', delimiter=',')
    #freqs_1000 = np.loadtxt(save_directory_1000 + '/freqs.txt', delimiter=',')
    #print(freqs_1000)


    alpha = 1
    plt.plot(distances_01a[::2], 2 * freqs_01a[::2, 0], color="green")
    plt.plot(distances_03[::2], 2 * freqs_03[::2, 0], color="gold")
    plt.plot(distances_05[::2], 2 * freqs_05[::2, 0], color="cyan")
    plt.plot(distances_07[::2], 2 * freqs_07[::2, 0], color="purple")


    plt.errorbar(distances_01a[::2], 2 * freqs_01a[::2, 0], freqs_01a[::2, 1], marker='o', ls='', ecolor="k", mec='black', color="green", label='$\\nu=0.020$', alpha=alpha)
    #plt.errorbar(distances_02, 2 * freqs_02[:, 0], freqs_02[:, 1], marker='o', ls='', ecolor="k", mec='black', color="purple", label='$\\nu=0.025$', alpha=alpha)
    plt.errorbar(distances_03[::2], 2 * freqs_03[::2, 0], freqs_03[::2, 1], marker='o', ls='', ecolor="k", mec='black', color="gold", label='$\\nu=0.030$', alpha=alpha)
    #plt.errorbar(distances_04, 2 * freqs_04[:, 0], freqs_04[:, 1], marker='o', ls='', ecolor="k", mec='black', color="red", label='$\\nu=0.035$', alpha=alpha)
    plt.errorbar(distances_05[::2], 2 * freqs_05[::2, 0], freqs_05[::2, 1], marker='o', ls='', ecolor="k", mec='black', color="cyan", label='$\\nu=0.040$', alpha=alpha)
    #plt.errorbar(distances_06, 2 * freqs_06[:, 0], freqs_06[:, 1], marker='o', ls='', ecolor="k", mec='black', color="blue", label='$\\nu=0.045$', alpha=alpha)
    plt.errorbar(distances_07[::2], 2 * freqs_07[::2, 0], freqs_07[::2, 1], marker='o', ls='', ecolor="k", mec='black', color="purple", label='$\\nu=0.050$', alpha=alpha)
    plt.hlines(0, 0, 60, colors="k")
    #plt.plot(distances_1000, freqs_1000[:, 0] / 200, color="purple", alpha=0.3)
    #plt.errorbar(distances_1000, freqs_1000[:, 0] / 200, freqs_1000[:, 1] / 200, marker='o', ls='', ecolor="k", mec='black', color="purple", label='$\\nu=0.32$', alpha=alpha)



    #plt.errorbar(distances_01a, 2 * freqs_01a[:, 0], freqs_01a[:, 1], marker='o', ls='', ecolor="k", mec='black', color="green", label='$\sigma_i=6\ \\textrm{mm}$', alpha=alpha)
    #plt.errorbar(distances_01b, 2 * freqs_01b[:, 0], freqs_01b[:, 1], marker='o', ls='', ecolor="k", mec='black', color="purple", label='$\sigma_i=9\ \\textrm{mm}$', alpha=alpha)
    #plt.errorbar(distances_01c, 2 * freqs_01c[:, 0], freqs_01c[:, 1], marker='o', ls='', ecolor="k", mec='black', color="gold", label='$\sigma_i=12\ \\textrm{mm}$', alpha=alpha)
    plt.xlabel('$d\ \\textrm{(mm)}$', fontsize=25)
    plt.ylabel('$\Omega \\times 10^{-3}$', fontsize=25)
    plt.xlim(20, 60)
    plt.ylim(-5, 100)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=15, loc="upper right")
    plt.grid(color="k", alpha=0.2)
    plt.tight_layout()
    plt.savefig("freqs.png", dpi=150)
    #plt.savefig("period_windows.png", dpi=200)
    plt.close()