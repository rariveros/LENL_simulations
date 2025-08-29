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

    # umbral
    thr1 = 0.01
    thr2 = 0.01

    # columnas 0 y 1 de modules
    m0 = powers[:, 0]
    m1 = powers[:, 1]


    # arranque en negro y enmascarar
    colors = np.full(m0.shape, "#5bf995", dtype=object)

    red_mask = (m0 < thr1) & (m1 < thr2)  # ambos < 0.001  -> rojo
    blue_mask = (m0 > thr1) & (m1 < thr2)  # ambos > 0.001  -> azul
    # lo demás queda negro

    colors[red_mask] = "#DF0A00"  # Damped
    colors[blue_mask] = "#0005f6"  # ROs

    freqs_masked = freqs.copy()
    freqs_masked[powers[:, 0] < 0.01] = 0

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3, 3.5))
    ax1.scatter(dists, 1000 * freqs_masked[:, 0], c=colors, alpha=1.0, edgecolors="k", lw=0.6, s=14, zorder=5)
    ax1.plot(dists, 1000 * freqs_masked[:, 0], color="k", lw=1)
    ax1.set_ylabel('$\Omega \\times 10^{-3}$', fontsize=15)
    ax1.tick_params(axis="both", direction="in", labeltop=False, labelbottom=False, top=False, bottom=False, labelsize=13)
    ax1.grid(alpha=0.3)
    ax1.hlines(0, 0, 40, colors="k")

    #ax2.scatter(dists, powers, s=70, marker="o", c="k", edgecolor="k", lw=0.8, zorder=5)
    #ax2.errorbar(dists, powers[:, 0], powers[:, 1], fmt="o", c="k")
    ax2.scatter(dists, powers[:, 0], c=colors, alpha=1.0, edgecolors="k", lw=0.6, s=14, zorder=5)
    ax2.plot(dists, powers[:, 0], color="k", lw=1)
    ax2.set_xlabel('$d$', fontsize=15)
    ax2.set_ylabel('$\\textrm{Power}$', fontsize=15)
    ax2.tick_params(axis="both", direction="in", labelsize=13)
    ax2.set_xticks([0, 10, 20, 30, 40])
    ax2.grid(alpha=0.3)
    ax2.hlines(0, 0, 40, colors="k")

    fig.subplots_adjust(left=0.23, right=0.95, bottom=0.27, top=0.95)
    plt.savefig('characterization.png', dpi=300)
    plt.close()