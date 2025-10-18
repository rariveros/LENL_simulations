from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    directory = r"C:\mnustes_science\simulation_data\FD\coherent\meanfield\Delta=0.0000\gamma=0.1000"
    Freqs = np.loadtxt(directory + '/frequencies.txt', delimiter=',')
    Ks = np.loadtxt(directory + '/ks.txt', delimiter=',')
    Omegas = np.loadtxt(directory + '/omegas.txt', delimiter=',')

    # plot 2D colormap
    # === Plot con mismo formato que parameter_space ===

    Freqs = filtro_superficie(Freqs, 2, "X")
    Freqs = np.nan_to_num(Freqs, nan=0.0)
    fig, ax = plt.subplots(1, 1, figsize=(6.0, 3.2))

    pc = ax.pcolormesh(Ks, Omegas, Freqs / 2, shading='auto', cmap='inferno', vmin=0, vmax=0.18)
    ax.tick_params(axis="both", direction="in", labelsize=11, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax.set_xlabel(r"$\kappa$", fontsize=15)
    ax.set_ylabel(r"$\Omega$", fontsize=15)
    #ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    #ax.set_yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25])

    # === Colorbar al costado, igual estilo ===
    cax = fig.add_axes([0.82, 0.15, 0.025, 0.75])
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label(r"$\omega$", fontsize=14)
    cb.ax.tick_params(labelsize=10)
    cb.set_ticks([0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18])

    # === Ajustes de márgenes ===
    fig.subplots_adjust(wspace=0.15, hspace=0.18, left=0.26, right=0.8, bottom=0.15, top=0.9)

    plt.savefig("freq_space_formatted.png", dpi=300)
    plt.close()