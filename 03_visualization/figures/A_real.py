from directories import *
from back_process import *


if __name__ == "__main__":
    disco = 'C'
    initial_dir_data = str(disco) + ':/Users/mnustes_science/PT_fluids/mnustes_science/simulation_data/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X_mm = np.loadtxt(directory + '/X.txt', delimiter=',')
    T_s = np.loadtxt(directory + '/T.txt', delimiter=',')

    t_i = int(3 * len(T_s) / 6)
    t_f = int(4 * len(T_s) / 6)
    #print(Z_mm[0, :])
    #plt.plot(X_mm, Z_mm[0, :])
    #plt.show()
    #plt.close()
    nu = 0.15
    #Z_mm = Z_mm[t_i:t_f, :]
    #T_s = T_s[t_i:t_f]
    Z_mm = np.zeros((len(T_s), len(X_mm)), dtype='complex_')
    for i in range(len(T_s)):
        for j in range(len(X_mm)):
            Z_mm[i, j] = (Z_r[i, j] + 1j * Z_i[i, j]) * np.exp(1j * np.sqrt(nu) * X_mm[j]) + (Z_r[i, j] - 1j * Z_i[i, j]) * np.exp(- 1 * 1j * np.sqrt(nu) * X_mm[j])
    np.savetxt(directory + '/A_real.txt', Z_mm, delimiter=',')

    fig, ax = plt.subplots()

    # legend
    pcm = plt.pcolormesh(X_mm, T_s, np.real(Z_mm), cmap=parula_map, vmin=np.amin(np.real(Z_mm)), vmax=np.amax(np.real(Z_mm)), shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$\psi(x, t)$', rotation=0, size=20, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlim([X_mm[0], X_mm[-1]])
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=20)

    #plt.ylim([20, 80])
    plt.ylabel('$t$', size='25')
    plt.yticks(fontsize=20)
    #plt.ylim([0, 0.5])

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/A_real.png', dpi=300)
    plt.close()