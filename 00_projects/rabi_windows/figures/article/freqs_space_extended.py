from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    disc = "X:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]

    FREQS = []
    K = []

    # --- solo cambio: dpi=300 (para impresión) ---
    # (si quieres tamaño controlado para 0.4 linewidth PRL: fig = plt.figure(figsize=(2.8, 2.2), dpi=300))
    fig = plt.figure(figsize=(5, 4), dpi=300)
    ax = fig.add_subplot(projection='3d')

    for directory_01 in directories:
        print(directory_01)
        dir_01 = working_directory + "/" + directory_01 + "/sigma=3.000/gamma=0.2800"
        freqs = np.loadtxt(dir_01 + '/freqs.txt', delimiter=',')
        powers = np.loadtxt(dir_01 + '/powers.txt', delimiter=',')
        Ks = np.loadtxt(dir_01 + '/dists.txt', delimiter=',')
        nus = np.loadtxt(dir_01 + '/nus.txt', delimiter=',')

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

        # tu scatter (edgecolor siempre negro)
        #ax.scatter(Ks[::2], nus[::2], 1000 * freqs_masked[::2, 0], c=colors[::2], alpha=1.0, edgecolors="k", lw=0.4, s=25, zorder=5)
        #ax.plot(Ks[::2], nus[::2], 1000 * freqs_masked[::2, 0], c="k", zorder=1)
        ax.scatter(Ks, nus, powers[:, 0], c=colors, alpha=1.0, edgecolors="k", lw=0.4, s=20, zorder=5)
        ax.plot(Ks, nus, powers[:, 0], c="k", zorder=1)
    # --- tu vista intacta ---
    #ax.view_init(23, -57, 0)
    ax.view_init(23, 60, 0)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis._axinfo['tick']['inward_factor'] = 0.0  # how far ticks go *into* the axes
        axis._axinfo['tick']['outward_factor'] = 0.3  # no outward ticks
    # ======== SOLO LO QUE PEDISTE: fondo blanco y sin grilla/paneles grises ========
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")
    ax.grid(True, which='major', color='0.5', linewidth=0.2, alpha=0.6)
    # paneles blancos, sin borde gris
    ax.xaxis.pane.set_facecolor("white")
    ax.yaxis.pane.set_facecolor("white")
    ax.zaxis.pane.set_facecolor("white")
    ax.xaxis.pane.set_edgecolor("white")
    ax.yaxis.pane.set_edgecolor("white")
    ax.zaxis.pane.set_edgecolor("white")
    # ==============================================================================

    # ======== Etiquetas y ticks legibles (activables) ========
    # (Descomenta si quieres etiquetas LaTeX y tamaños consistentes para figura pequeña)
    ax.set_xlabel(r'$d$', fontsize=20, labelpad=8)   # EJE K = kappa
    ax.set_ylabel(r'$\nu$', fontsize=20, labelpad=8)      # EJE NU = nu
    ax.set_zlabel(r'$\textrm{Power}$', fontsize=20, labelpad=8)   # FRECUENCIA = Omega

    # Ticks más legibles sin cambiar tu data
    ax.tick_params(axis='both', which='both', direction="in", labelsize=13, length=6, width=1.6, pad=1)
    ax.tick_params(axis='z',    which='both', direction="in", labelsize=13, length=6, width=1.6, pad=1)
    # ========================================================

    # --- tus límites y ticks están acá por si quieres reactivarlos ---
    #ax.set_xlim(-0.5, 8)
    ax.set_ylim(0.23, 0.33)
    ax.set_zlim(0, 3)
    # ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], ["$0.0$", "$0.1$", "$0.2$", "$0.3$", "$0.4$", "$0.5$"], fontsize=15)
    # ax.set_yticks([0.0, 0.05, 0.1, 0.150, 0.2], ["$0.00$", "$0.05$", "$0.10$", "$0.15$", "$0.20$"], fontsize=15)
    ax.set_zticks([0, 1, 2, 3])
    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.05, top=0.98)
    ax.set_box_aspect([2, 2, 0.8])
    # --- guardado (si quieres) ---
    plt.savefig("freq_space.png", dpi=300)
    # plt.savefig("freq_space_01.pdf", bbox_inches='tight')

    #plt.show()