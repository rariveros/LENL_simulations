from back_process import *


if __name__ == "__main__":
    disco = 'C:/'
    initial_dir_data = str(disco) + ':/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z = np.loadtxt('C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/drift_fernando/OP_Mod_a90.dat')
    X = (np.arange(len(Z[0, :])) - len(Z[0, :]) / 2) * 0.25
    T = np.arange(len(Z[:, 0])) * 0.006

    fig, ax = plt.subplots()
    file = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/drift_fernando"
    subfile = "/OP_Mod_a90"
    if not os.path.exists(file + subfile):
        os.makedirs(file + subfile)
    np.savetxt(file + subfile + '/Z.txt', Z, delimiter=',')
    np.savetxt(file + subfile + '/X.txt', X, delimiter=',')
    np.savetxt(file + subfile + '/T.txt', T, delimiter=',')

    pcm = plt.pcolormesh(X, T, Z, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$|A|$', rotation=0, size=25, labelpad=-27, y=1.11)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    plt.tight_layout()
    plt.savefig(file + subfile + '/field_visual.png', dpi=200)
    plt.close()