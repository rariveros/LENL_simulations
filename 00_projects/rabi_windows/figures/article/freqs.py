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
    freqs = np.loadtxt(directory + '/analysis/freqs.txt', delimiter=',')
    dists = np.loadtxt(directory + '/analysis/dists.txt', delimiter=',')
    print(dists)
    print(freqs)
    plt.scatter(dists, freqs[:, 0], c="k")
    plt.xlabel("$d$", fontsize=18)
    plt.ylabel("$\Omega$", fontsize=18)
    plt.show()