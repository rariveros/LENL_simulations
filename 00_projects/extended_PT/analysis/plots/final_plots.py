from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    #data = [parameter, period, period_std, power_L, power_R]
    data = np.loadtxt(directory + '/data.txt', delimiter=',')
    data_dist = np.loadtxt(directory + '/data_dist.txt', delimiter=',')
    #plt.scatter(data[1:-2, 0], 100/data[1:-2, 1], c="k", s=45)
    print(data[1:-2, 1])
    plt.errorbar(data[1:-3, 0], 100 / data[1:-3, 1], 100 * np.abs(data[1:-3, 2]/data[1:-3, 1]**2), marker='o', ls='',
                 ecolor="k", mec='black',
                 color="k")
    plt.xlabel('$\gamma_0$', size='25')
    plt.xticks(fontsize="20")
    plt.ylabel('$f_{\\textrm{RO}}$', size='25')
    plt.yticks(fontsize="20")
    plt.grid(linestyle='--', alpha=0.5)
    plt.title("$\\times 10^{-2}$", loc='left', size="18")
    plt.tight_layout()
    plt.savefig(directory + '/figures/freq.png', dpi=300)
    plt.close()

    dist = np.append(data_dist[1:-3, 0], [29.505, 30, 31, 32, 33, 34, 35])
    freq_dist = np.append(100 / data_dist[1:-3, 1], [0, 0, 0, 0, 0, 0, 0])
    freq_dist_err = np.append(100 * np.abs(data_dist[1:-3, 2]/data_dist [1:-3, 1]**2), [0, 0, 0, 0, 0, 0, 0])
    plt.fill_between([22, 29.505], -0.1, 1, color=(194 / 255, 226 / 255, 191 / 255), alpha=0.8, zorder=0)
    plt.fill_between([29.505, 35], -0.1, 1, color=(191 / 255, 192 / 255, 225 / 255), alpha=0.8, zorder=0)

    plt.errorbar(dist, freq_dist, freq_dist_err, marker='o', ls='', ecolor="k", mec='black', color="k")
    plt.xlabel('$d\ \\textrm{(mm)}$', size='25')
    plt.xticks(fontsize="20")
    plt.ylabel('$f_{\\textrm{RO}}$', size='25')
    plt.yticks(fontsize="20")
    plt.ylim([0, 1])
    plt.xlim([23, 35])
    plt.text(23.4, 0.22, "$\\textrm{Rabi Oscillations}$", fontsize=20)
    plt.text(30.2, 0.5, "$\\textrm{Non-interacting}$", fontsize=20)
    plt.grid(linestyle='--', alpha=0.5)
    plt.title("$\\times 10^{-2}$", loc='left', size="18")
    plt.tight_layout()
    plt.savefig(directory + '/figures/freq_dist.png', dpi=300)
    plt.show()
    plt.close()


    plt.scatter(data[1:-3, 0], data[1:-3, 3], c="k", s=45)
    plt.xlabel('$\gamma_0$', size='25')
    plt.xticks(fontsize="20")
    plt.ylabel('$P_{\\textrm{RO}}$', size='25')
    plt.yticks(fontsize="20")
    plt.grid(linestyle='--', alpha=0.5)
    plt.title("$\\times 10^{-2}$", loc='left', size="18", color="white")
    plt.tight_layout()
    plt.savefig(directory + '/figures/power.png', dpi=300)
    plt.close()