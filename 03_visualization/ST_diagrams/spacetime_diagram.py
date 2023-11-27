from back_process import *


if __name__ == "__main__":
    disco = 'C:/'
    initial_dir_data = str(disco) + ':/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')

    beta = 0.004811649356064012
    ####

    Z_complex = Z_r + 1j * Z_i
    arg = np.angle(Z_complex)

    # Definiendo variables finales
    Z_modulo = np.absolute(Z_complex)
    Z_modulo = filtro_superficie(Z_modulo, 5, "XY")
    fig, ax1 = plt.subplots()
    ti, tf = T[-401], T[-1]
    xi, xf = -100, 100

    #matplotlib.rc('xtick', labelsize=20)
    #matplotlib.rc('ytick', labelsize=20)

    #plt.plot(X, arg[ti, :])
    #plt.show()

    # legend
    pcm = ax1.pcolormesh(X, T - ti, np.flip(Z_modulo / np.sqrt(beta), axis=1), cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    #cbar.set_label('$|A(x,t)|$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell

    ax1.set_xlim([xi, xf])
    ax1.set_xlabel('$x$', size='25')
    ax1.tick_params(labelsize=18)

    ax1.set_ylim([0, tf - ti])
    ax1.set_ylabel('$t$', size='25')
    #plt.set_yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/field_visual.png', dpi=300)
    plt.close()