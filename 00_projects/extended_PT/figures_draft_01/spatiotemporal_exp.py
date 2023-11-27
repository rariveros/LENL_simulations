from directories import *
from back_process import *


if __name__ == "__main__":
    disco = 'N:'
    initial_dir_data = str(disco) + ':/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z = np.loadtxt(directory + '/Z_strobo.txt', delimiter=',')
    X = np.loadtxt(directory + '/X_mm.txt', delimiter=',')
    T = np.loadtxt(directory + '/T_strobo.txt', delimiter=',')

    Z = filtro_superficie(Z, 5, "YX")

    file_name = os.path.basename(directory)
    name_list = file_name.split("_")
    notation = 'fa'
    if notation == 'fa':
        f_i = float(name_list[0].split("=")[-1])
        a = float(name_list[1].split("=")[-1])
    elif notation == 'af':
        a = float(name_list[0].split("=")[-1])
        f_i = float(name_list[1].split("=")[-1])
    alpha, beta, nu, gamma = fluid_pdnls_parameters(f_i, a, d=20)

    # FIXES
    #X = X[8:-9]
    fig, ax = plt.subplots()
    ti, tf = 0, 400
    xi, xf = -150, 150
    x0 = -9

    g = 9790
    l_y = 16
    d = 20
    k_y = np.pi / l_y
    k = k_y
    tau = np.tanh(k * d)
    w_1 = np.sqrt(g * k * tau)

    # legend
    fig, ax1 = plt.subplots()
    pcm = ax1.pcolormesh(X - x0, T * (w_1 / (2 * np.pi)), Z, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    #cbar.set_label('$|A(x,t)|$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell

    ax1.set_xlim([xi, xf])
    ax1.set_xlabel('$x$', size='25')
    ax1.tick_params(labelsize=18)

    ax1.set_ylim([ti, tf])
    ax1.set_ylabel('$t$', size='25')
    # plt.set_yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/field_visual.png', dpi=300)
    plt.close()