from back_process import *

if __name__ == "__main__":
    disco = 'C:/'
    initial_dir_data = str(disco) + ':/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    #Z = np.loadtxt('C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/drift_fernando/OP_Mod_a90.dat')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    params = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    #[alpha, beta, gamma_0, dist,  mu, nu, sigma]
    dist = params[3]
    beta = 0.004811649356064012

    g = 9790
    l_y = 16
    d = 20
    k_y = np.pi / l_y
    k = k_y
    tau = np.tanh(k * d)
    w_1 = np.sqrt(g * k * tau)

    T = T / (w_1 / np.pi)

    Z_complex = Z_r + 1j * Z_i
    arg = np.angle(Z_complex)

    # Definiendo variables finales
    Z = np.absolute(Z_complex)

    fig, ax = plt.subplots()
    ti, tf = 300, 4999
    xi, xf = -100, 100

    # legend
    pcm = plt.pcolormesh(X, T - ti, Z / np.sqrt(beta), cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$|A|$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    #plt.xlim([xi, xf])
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylim([0, 54])
    plt.ylabel('$t\ \\textrm{(s)}$', size='25')
    plt.yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(str(dist) + '.png', dpi=300)
    plt.close()