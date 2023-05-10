from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    freq = np.loadtxt(directory + '/freq.txt', delimiter=',')
    ampd = np.loadtxt(directory + '/ampd.txt', delimiter=',')
    param = np.loadtxt(directory + '/param.txt', delimiter=',')

    plt.scatter(param, freq, c="k", zorder=10)
    plt.grid(alpha=0.5, zorder=1)
    plt.ylim([0, 1.1 * np.amax(freq)])
    plt.xlabel("$\gamma_0$", fontsize=20)
    plt.ylabel("$f$", fontsize=20)
    plt.tight_layout()
    plt.savefig(directory + '/freq_oscillations.png', dpi=300)
    plt.close()

    plt.scatter(param, ampd, c="k", zorder=10)
    plt.grid(alpha=0.5, zorder=1)
    plt.ylim([0, 1.1 * np.amax(freq)])
    plt.xlabel("$\gamma_0$", fontsize=20)
    plt.ylabel("$\delta A_R$", fontsize=20)
    plt.tight_layout()
    plt.savefig(directory + '/ampd_oscillations.png', dpi=300)
    plt.close()