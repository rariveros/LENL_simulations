import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    disco = 'C:/'
    initial_dir_data = str(disco) + ':/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    X_max = np.loadtxt(directory + '/X_max_sigma.txt', delimiter=',')
    sigmas = np.arange(4, 25)

    fig, ax = plt.subplots()
    textstr = '\n'.join((
        "$\\nu=-0.3$",
        "$\gamma_0 = 0.15$",
        "$\mu = 0.1$"))
    plt.scatter(sigmas, X_max, s=45, c="k", zorder=10)
    plt.grid(alpha=0.2, zorder=1)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("$\sigma_i\ (\\textrm{mm})$", fontsize=25)
    plt.ylabel("$x_{\\textrm{s}}\ (\\textrm{mm})$", fontsize=25)
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.77, 0.95, textstr, transform=ax.transAxes, fontsize=18,
            verticalalignment='top', bbox=props)
    plt.tight_layout()
    plt.savefig("Xmax_sigma.png", dpi=200)