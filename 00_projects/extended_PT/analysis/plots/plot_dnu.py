from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    parameters = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    index_01 = np.loadtxt(directory + '/index_01.txt', delimiter=',')
    index_02 = np.loadtxt(directory + '/index_02.txt', delimiter=',')
    index_03_real = np.loadtxt(directory + '/index_03_real.txt', delimiter=',')
    index_03_imag = np.loadtxt(directory + '/index_03_imag.txt', delimiter=',')
    freqs = np.loadtxt(directory + '/freqs.txt', delimiter=',')
    phases = np.loadtxt(directory + '/phases.txt', delimiter=',')

    x_left = np.arange(0, 0.081, 0.001)
    x_center = np.arange(0.08, 0.17, 0.001)
    x_right = np.arange(0.17, 0.211, 0.001)

    fig, ax = plt.subplots()
    plt.fill_between(x_left, 0, 1.1 * np.amax(index_01), color=(14 / 255, 90 / 255, 111 / 255), alpha=0.2, zorder=0)
    plt.vlines(0.08, 0, 1.1 * np.amax(index_01), colors="k", linestyle="--")
    plt.fill_between(x_center, 0, 1.1 * np.amax(index_01), color=(70 / 255, 169 / 255, 64 / 255), alpha=0.2, zorder=0)
    plt.vlines(0.17, 0, 1.1 * np.amax(index_01), colors="k", linestyle="--")
    plt.fill_between(x_right, 0, 1.1 * np.amax(index_01), color=(78 / 255, 0 / 255, 136 / 255), alpha=0.2, zorder=0)
    plt.scatter(parameters, index_01, c="k")
    plt.xlabel('$\\nu_i$', size='20')
    plt.ylabel('$\\xi^{1}_s$', size='20')
    plt.xlim(0, 0.21)
    plt.ylim(0, 1.1 * np.amax(index_01))
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(directory + '//index_01.png', dpi=300)
    plt.close()

    fig, ax = plt.subplots()
    plt.fill_between(x_left, 0, 1.1 * np.amax(index_02), color=(14 / 255, 90 / 255, 111 / 255), alpha=0.2, zorder=0)
    plt.vlines(0.08, 0, 1.1 * np.amax(index_02), colors="k", linestyle="--")
    plt.fill_between(x_center, 0, 1.1 * np.amax(index_02), color=(70 / 255, 169 / 255, 64 / 255), alpha=0.2, zorder=0)
    plt.vlines(0.17, 0, 1.1 * np.amax(index_02), colors="k", linestyle="--")
    plt.fill_between(x_right, 0, 1.1 * np.amax(index_02), color=(78 / 255, 0 / 255, 136 / 255), alpha=0.2, zorder=0)
    plt.scatter(parameters, index_02, c="k")
    plt.xlabel('$\\nu_i$', size='20')
    plt.ylabel('$\\xi^{2}_s$', size='20')
    plt.xlim(0, 0.21)
    plt.ylim(0, 1.1 * np.amax(index_02))
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(directory + '//index_02.png', dpi=300)
    plt.close()

    fig, ax = plt.subplots()
    plt.fill_between(x_left, 1.1 * np.amin(index_03_imag), 1.1 * np.amax(index_03_real), color=(14 / 255, 90 / 255, 111 / 255), alpha=0.2, zorder=0)
    plt.vlines(0.08, 1.1 * np.amin(index_03_imag), 1.1 * np.amax(index_03_real), colors="k", linestyle="--")
    plt.fill_between(x_center, 1.1 * np.amin(index_03_imag), 1.1 * np.amax(index_03_real), color=(70 / 255, 169 / 255, 64 / 255), alpha=0.2, zorder=0)
    plt.vlines(0.17, 1.1 * np.amin(index_03_imag), 1.1 * np.amax(index_03_real), colors="k", linestyle="--")
    plt.fill_between(x_right, 1.1 * np.amin(index_03_imag), 1.1 * np.amax(index_03_real), color=(78 / 255, 0 / 255, 136 / 255), alpha=0.2, zorder=0)
    plt.scatter(parameters, index_03_real, c="r", edgecolors='black', label="$\\textrm{Re}\ (\\zeta_{\\textrm{sym}})$", zorder=10)
    plt.scatter(parameters, index_03_imag, c="b", edgecolors='black', label="$\\textrm{Im}\ (\\zeta_{\\textrm{sym}})$", zorder=10)
    plt.xlabel('$\\nu$', size='25')
    plt.ylabel('$\\zeta_{\\textrm{sym}}$', size='25')
    plt.hlines([0, 1], 1.1 * np.amin(index_03_imag), 1.1 * np.amax(index_03_real), colors=["k", "k"], zorder=2)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim(0, 0.2)
    plt.ylim(1.1 * np.amin(index_03_imag), 1.1 * np.amax(index_03_real))
    plt.grid(linestyle='--', alpha=0.5, zorder=1)
    plt.legend(loc="center left", fontsize=15)
    plt.tight_layout()
    plt.savefig(directory + '//index_03.png', dpi=300)
    plt.close()

    fig, ax = plt.subplots()
    plt.errorbar(parameters, (freqs[:, 0] + freqs[:, 2]) / 2, (freqs[:, 1] + freqs[:, 3]) / 2, marker='o', ls='', ecolor="k", mec='black',
                 color="k")
    plt.fill_between(x_left, 0, 0.031, color=(14 / 255, 90 / 255, 111 / 255), alpha=0.2, zorder=0)
    plt.vlines(0.08, 0, 0.031, colors="k", linestyle="--")
    plt.fill_between(x_center, 0, 0.031, color=(70 / 255, 169 / 255, 64 / 255), alpha=0.2, zorder=0)
    plt.vlines(0.17, 0, 0.031, colors="k", linestyle="--")
    plt.fill_between(x_right, 0, 0.031, color=(78 / 255, 0 / 255, 136 / 255), alpha=0.2, zorder=0)
    plt.xlabel('$\\nu$', size='25')
    plt.ylabel('$f_{\\textrm{so}}$', size='25')
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim(0, 0.21)
    plt.ylim(0.0, 0.03)
    plt.grid(linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(directory + '//freqs.png', dpi=300)
    plt.close()

    fig, ax = plt.subplots()
    ax.errorbar(parameters, phases[:, 0], phases[:, 1], marker='o', ls='', ecolor="k", mec='black',
                 color=(216 / 255, 135 / 255, 25 / 255), label="$x_L$")
    ax.errorbar(parameters, phases[:, 2], phases[:, 3], marker='o', ls='', ecolor="k", mec='black',
                 color=(49 / 255, 0 / 255, 129 / 255), label="$x_R$")
    plt.fill_between(x_left, 0, 2, color=(14 / 255, 90 / 255, 111 / 255), alpha=0.2, zorder=0)
    plt.vlines(0.08, 0, 2, colors="k", linestyle="--")
    plt.fill_between(x_center, 0, 2, color=(70 / 255, 169 / 255, 64 / 255), alpha=0.2, zorder=0)
    plt.vlines(0.17, 0, 2, colors="k", linestyle="--")
    plt.fill_between(x_right, 0, 2, color=(78 / 255, 0 / 255, 136 / 255), alpha=0.2, zorder=0)
    plt.xlabel('$\\nu$', size='20')
    plt.ylabel('$\phi$', size='20', rotation=0)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim(0, 0.21)
    plt.ylim(0, 2)
    ax.yaxis.set_major_formatter(FormatStrFormatter('$%g \pi$'))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=0.25))
    plt.grid(linestyle='--', alpha=0.5)
    plt.legend(loc="upper right", fontsize=20)
    plt.savefig(directory + '//phases.png', dpi=300)
    plt.close()
