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


    x_right = np.arange(9.5, 26, 0.01)

    fig, ax = plt.subplots()
    #plt.fill_between(x_right, 0, 1.5, color=(78 / 255, 0 / 255, 136 / 255), alpha=0.2, zorder=0)
    plt.scatter(parameters, index_03_real, c="r", edgecolors='black', label="$\\textrm{Re}\ \\xi_S$")
    plt.scatter(parameters, index_03_imag, c="b", edgecolors='black', label="$\\textrm{Im}\ \\xi_S$")
    plt.xlabel('$a$', size='20')
    plt.ylabel('$\\xi_S$', size='20')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #plt.xlim(13.75, 24.5)
    #plt.ylim(0, 1.4)
    plt.grid(linestyle='--', alpha=0.5)
    plt.legend(loc="upper right", fontsize=20)
    plt.tight_layout()
    plt.savefig(directory + '//index_03.png', dpi=300)
    plt.close()

    fig, ax = plt.subplots()
    plt.errorbar(parameters, freqs[:, 0], freqs[:, 1], marker='o', ls='', ecolor="k", mec='black',
                 color=(216 / 255, 135 / 255, 25 / 255), label="$x_L$")
    plt.errorbar(parameters, freqs[:, 2], freqs[:, 3], marker='o', ls='', ecolor="k", mec='black',
                 color=(49 / 255, 0 / 255, 129 / 255), label="$x_R$")
    #plt.fill_between(x_right, 0.0, 0.04, color=(78 / 255, 0 / 255, 136 / 255), alpha=0.2, zorder=0)
    plt.xlabel('$a$', size='20')
    plt.ylabel('$\\frac{\omega}{2 \pi}$', size='20', rotation=0)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #plt.xlim(13.75, 24.5)
    #plt.ylim(0.0, 0.04)
    plt.grid(linestyle='--', alpha=0.5)
    plt.legend(loc="upper right", fontsize=20)
    plt.tight_layout()
    plt.savefig(directory + '//freqs.png', dpi=300)
    plt.close()

    fig, ax = plt.subplots()
    ax.errorbar(parameters, phases[:, 0], phases[:, 1], marker='o', ls='', ecolor="k", mec='black',
                 color=(216 / 255, 135 / 255, 25 / 255), label="$x_L$")
    ax.errorbar(parameters, phases[:, 2], phases[:, 3], marker='o', ls='', ecolor="k", mec='black',
                 color=(49 / 255, 0 / 255, 129 / 255), label="$x_R$")
    #plt.fill_between(x_right, 0, 2, color=(78 / 255, 0 / 255, 136 / 255), alpha=0.2, zorder=0)
    plt.xlabel('$a$', size='20')
    plt.ylabel('$\phi$', size='20', rotation=0)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #plt.xlim(13.75, 24.5)
    #plt.ylim(0, 2)
    ax.yaxis.set_major_formatter(FormatStrFormatter('$%g \pi$'))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=0.25))
    plt.grid(linestyle='--', alpha=0.5)
    plt.legend(loc="upper right", fontsize=20)
    plt.savefig(directory + '//phases.png', dpi=300)
    plt.close()

