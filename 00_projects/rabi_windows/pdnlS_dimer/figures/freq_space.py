from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    disc = "D:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]

    FREQS = []
    K = []

    # --- solo cambio: dpi=300 (para impresión) ---
    # (si quieres tamaño controlado para 0.4 linewidth PRL: fig = plt.figure(figsize=(2.8, 2.2), dpi=300))
    fig = plt.figure(figsize=(6, 4.5), dpi=300)
    ax = fig.add_subplot(projection='3d')

    for directory_01 in directories:
        print(directory_01)
        dir_01 = working_directory + "/" + directory_01 + "/mu=0.1000/gamma=0.2000"
        freqs = np.loadtxt(dir_01 + '/frequencies.txt', delimiter=',')
        modules = np.loadtxt(dir_01 + '/modules.txt', delimiter=',')
        Ks = np.loadtxt(dir_01 + '/Ks.txt', delimiter=',')
        nus = np.loadtxt(dir_01 + '/nus.txt', delimiter=',')

        # umbral
        # umbral
        thr1 = 0.01
        thr2 = 0.01

        # columnas 0 y 1 de modules
        m0 = modules[:, 0]
        m1 = modules[:, 1]

        # arranque en negro y enmascarar
        colors = np.full(m0.shape, "#5bf995", dtype=object)

        red_mask = (m0 < thr1) & (m1 < thr2)  # ambos < 0.001  -> rojo
        blue_mask = (m0 > thr1) & (m1 < thr2)  # ambos > 0.001  -> azul
        # lo demás queda negro

        colors[red_mask] = "#DF0A00"  # Damped
        colors[blue_mask] = "#0005f6"  # ROs
        # tu scatter (edgecolor siempre negro)
        ax.scatter(100 * Ks, nus, 1000 * freqs, c=colors, zorder=4, alpha=1.0, edgecolors="k", lw=0.6, s=20)
        ax.plot(100 * Ks, nus, 1000 * freqs, c="k", zorder=1)
    # --- tu vista intacta ---
    ax.view_init(30, -167, 0)

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
    ax.set_xlabel(r'$\kappa \times 10^{-2}$', fontsize=20)   # EJE K = kappa
    ax.set_ylabel(r'$\nu$', fontsize=20)      # EJE NU = nu
    ax.set_zlabel(r'$\Omega \times 10^{-3}$', fontsize=20)   # FRECUENCIA = Omega

    # Ticks más legibles sin cambiar tu data
    ax.tick_params(axis='both', which='major', labelsize=15, length=6, width=1.6, pad=2)
    ax.tick_params(axis='z',    which='major', labelsize=15, length=6, width=1.6, pad=2)
    # ========================================================

    # --- tus límites y ticks están acá por si quieres reactivarlos ---
    ax.set_xlim(-0.5, 8)
    ax.set_ylim(-0.01, 0.265)
    ax.set_zlim(0, 30)
    # ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], ["$0.0$", "$0.1$", "$0.2$", "$0.3$", "$0.4$", "$0.5$"], fontsize=15)
    # ax.set_yticks([0.0, 0.05, 0.1, 0.150, 0.2], ["$0.00$", "$0.05$", "$0.10$", "$0.15$", "$0.20$"], fontsize=15)
    # ax.set_zticks([0.0, 0.05, 0.1, 0.150], ["$0.00$", "$0.05$", "$0.10$", "$0.15$"], fontsize=15)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.2, top=0.98)
    ax.set_box_aspect([2, 2, 0.8])
    # --- guardado (si quieres) ---
    plt.savefig("freq_space.png", dpi=300)
    # plt.savefig("freq_space_01.pdf", bbox_inches='tight')

    #plt.show()