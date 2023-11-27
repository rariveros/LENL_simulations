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

    xmin, xmax = 0.1875, 0.285
    ymin, ymax = 1.26, 1.52
    fig, ax = plt.subplots()
    plt.errorbar(data[1:-3, 0], 100 / data[1:-3, 1], 100 * np.abs(data[1:-3, 2]/data[1:-3, 1]**2), marker='o', ls='',
                 ecolor="k", mec='black',
                 color="k")
    #ax.errorbar(data[1:-3, 0], data[1:-3, 1], data[1:-3, 2], marker='o', ls='',
    #             ecolor="k", mec='black',
    #             color="k")
    ax.fill_between([0.17, 0.225], ymin, ymax, color=(194 / 255, 226 / 255, 191 / 255), alpha=0.8, zorder=0)
    ax.fill_between([0.2251, 0.3], ymin, ymax, color=(240 / 255, 215 / 255, 180/ 255), alpha=0.8, zorder=0)
    ax.vlines(0.225, ymin, ymax, colors="k", linewidth=1, linestyles="--")
    ax.set_aspect(0.42 * (xmax - xmin) / (ymax - ymin))

    plt.xlabel('$\gamma_0$', size='25')
    plt.xticks(fontsize="20")
    plt.ylabel('$f_{\\textrm{RO}}$', size='25')
    plt.yticks(fontsize="20")
    plt.grid(linestyle='--', alpha=0.5)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.title("$\\times 10^{-2}$", loc='left', size="18")
    plt.tight_layout()
    plt.savefig(directory + '/figures/freq.png', dpi=300)
    plt.close()

    xmin, xmax = 23, 38
    ymin, ymax = -0.05, 1
    fig, ax = plt.subplots()
    dist = np.append(data_dist[:, 0], np.arange(30.5, 38, 0.25))
    freq_dist = np.append(100 / data_dist[:, 1], 0*np.arange(30.5, 38, 0.25))
    freq_dist_err = np.append(100 * np.abs(data_dist[:, 2]/data_dist[:, 1]**2), 0*np.arange(30.5, 38, 0.25))
    ax.fill_between([22, 30.5], ymin, ymax, color=(194 / 255, 226 / 255, 191 / 255), alpha=0.8, zorder=0)
    ax.fill_between([30.5, 40], ymin, ymax, color=(191 / 255, 192 / 255, 225 / 255), alpha=0.8, zorder=0)
    ax.hlines(0, 20, 41, colors="k", linewidth=1)
    ax.vlines(30.5, -0.1, 1, colors="k", linewidth=1, linestyles="--")
    ax.errorbar(dist, freq_dist, freq_dist_err, marker='o', ls='', ecolor="k", mec='black', color="k")
    ax.set_aspect(1 * (xmax - xmin) / (ymax - ymin))

    plt.xlabel('$d\ \\textrm{(mm)}$', size='20')
    plt.xticks(fontsize="16")
    plt.ylabel('$f_{\\textrm{RO}}$', size='20')
    plt.yticks(fontsize="16")
    plt.ylim([ymin, ymax])
    plt.xlim([xmin, xmax])
    plt.text(23.2, 0.1, "$\\textrm{Rabi Oscillations}$", fontsize=16)
    plt.text(31, 0.1, "$\\textrm{Non-interacting}$", fontsize=16)
    plt.grid(linestyle='--', alpha=0.5)
    plt.title("$\\times 10^{-2}$", loc='left', size="18")
    plt.tight_layout()
    plt.savefig(directory + '/figures/freq_dist.png', dpi=350)
    #plt.show()
    plt.close()


    plt.scatter(data[1:-3, 0], data[1:-3, 3], c="k", s=45)
    plt.xlabel('$\gamma_0$', size='25')
    plt.xticks(fontsize="20")
    plt.ylabel('$P_{\\textrm{RO}}$', size='25')
    plt.yticks(fontsize="20")
    plt.grid(linestyle='--', alpha=0.5)
    plt.title("$\\times 10^{-2}$", loc='left', size="18", color="white")
    plt.tight_layout()
    plt.savefig(directory + '/figures/power.png', dpi=190)
    plt.close()