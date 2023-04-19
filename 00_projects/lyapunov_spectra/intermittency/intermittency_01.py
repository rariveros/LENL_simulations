from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + "mnustes_science/simulation_data/FD/PDNLS_chaosa"
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    H = np.loadtxt(directory + '/Hamiltonian_t.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')[0:-1]

    H_filt = filtro_array(300, H)
    plt.plot(T, H, zorder=0, alpha=0.5, color="r", label="$H(t)$")
    plt.plot(T, H_filt, zorder=1, color="k", label="$\langle H(t) \\rangle_{t}$")
    plt.legend()
    plt.xlabel('$t$', size='20')
    plt.xticks(fontsize=15)
    plt.xlim([T[0], T[-1]])
    plt.ylabel('$H(t)$', size='20')
    plt.yticks(fontsize=15)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(directory + '/hamiltonian_mean.png', dpi=300)
    plt.close()