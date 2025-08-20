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
    ### E:\mnustes_science\simulation_data\FD\rabi_windows\dimensionless\var_dist\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800

    freqs = np.loadtxt(directory + '/freqs.txt', delimiter=',')
    powers = np.loadtxt(directory + '/powers.txt', delimiter=',')
    dists = np.loadtxt(directory + '/dists.txt', delimiter=',')

    # Make a copy if you don't want to overwrite original freqs
    freqs_masked = freqs.copy()
    print(frequencies)
    print(dists)
    # Apply condition: set freqs to zero where powers < 0.02
    freqs_masked[powers < 0.02] = 0
    print(len(dists))
    print(len(powers))
    fig, ax1 = plt.subplots()
    #ax1.scatter(dists, 1000 * freqs_masked[:, 0], s=70, marker="o", c="k", edgecolor="k", lw=0.8, zorder=5)
    ax1.scatter(dists, powers, s=70, marker="o", c="k", edgecolor="k", lw=0.8, zorder=5)
    ax1.set_xlabel('$d\ \\textrm{(mm)}$', fontsize=25)
    ax1.set_ylabel('$\Omega\ \\textrm{(mHz)}$', fontsize=25)
    ax1.tick_params(labelsize=20)
    plt.show()