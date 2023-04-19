from back_process import *


if __name__ == "__main__":
    disco = 'C:/'
    initial_dir_data = str(disco) + ':/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    Nt = len(T)

    ####

    Z_complex = Z_r + 1j * Z_i

    # Definiendo variables finales
    Z_modulo= np.absolute(Z_complex)

    fig, ax = plt.subplots()
    ti, tf = 0, 1000 #int(Nt * 0.5), int(Nt * (350 / 500))
    xi, xf = -100, 100

    # legend
    pcm = plt.pcolormesh(X, T, Z_modulo, cmap=parula_map, vmin=np.amin(Z_modulo), vmax=np.amax(Z_modulo), shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$|A|$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlim([xi, xf])
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylim([ti, tf])
    plt.ylabel('$t$', size='25')
    plt.yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/field_visual.png', dpi=300)
    plt.close()