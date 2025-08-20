from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    frequencies = []
    disc = 'D:/'

    initial_dir_data = str(disc) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    ### D:\mnustes_science\simulation_data\FD\PT_dimer\alpha=1.000\beta=1.000\mu=0.100\nu=0.320\sigma=3.000\gamma=0.280

    freqs = np.loadtxt(directory + '/analysis/freqs.txt', delimiter=',')
    powers = np.loadtxt(directory + '/analysis/powers.txt', delimiter=',')
    modules = np.loadtxt(directory + '/analysis/modules.txt', delimiter=',')
    dists = np.loadtxt(directory + '/analysis/dists.txt', delimiter=',')

    # Make a copy if you don't want to overwrite original freqs
    freqs_masked = freqs.copy()
    print(frequencies)
    print(dists)
    # Apply condition: set freqs to zero where powers < 0.02
    freqs_masked[powers < 0.02] = 0
    print(len(dists))
    print(len(powers))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 4), dpi=100)
    ax1.scatter(dists, 1000 * freqs_masked[:, 0], s=70, marker="o", c="k", edgecolor="k", lw=0.8, zorder=5)
    ax1.set_ylabel('$\Omega\ \\textrm{(mHz)}$', fontsize=25)
    ax1.tick_params(axis="both", direction="in", labeltop=False, labelbottom=False, top=False, bottom=False, labelsize=20)

    #ax2.scatter(dists, powers, s=70, marker="o", c="k", edgecolor="k", lw=0.8, zorder=5)
    ax2.errorbar(dists, modules[:, 0], modules[:, 1], fmt="o", c="k")
    ax2.errorbar(dists, modules[:, 2], modules[:, 3], fmt="o", c="k")
    ax2.set_xlabel('$d\ \\textrm{(mm)}$', fontsize=25)
    ax2.set_ylabel('$\\textrm{Power}$', fontsize=25)
    ax2.tick_params(axis="both", direction="in", labelsize=20)
    plt.savefig(directory + '/analysis/characterization.png', dpi=300)
    plt.show()
    plt.close()