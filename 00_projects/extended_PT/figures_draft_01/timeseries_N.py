import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FDa'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    print("#############   " + directory + "   #############")
    Z_modulo = np.loadtxt(directory + '/Z_strobo.txt', delimiter=',')
    X = np.loadtxt(directory + '/X_mm.txt', delimiter=',')
    T = np.loadtxt(directory + '/T_strobo.txt', delimiter=',')
    injection = np.loadtxt(directory + '/IL_mm.txt', delimiter=',')


    Nx = len(X)
    Nt = len(T)

    g = 9790
    l_y = 16
    d = 20
    k_y = np.pi / l_y
    k = k_y
    tau = np.tanh(k * d)
    w_1 = np.sqrt(g * k * tau)

    T = T * (w_1 / (2 * np.pi))

    L = 200
    x0 = -2
    #x0 = -9
    #x0 = 0

    inj_L = injection[1] - x0
    inj_R = injection[0] - x0

    d = "dist_01"
    middle = np.argmin(np.abs(X - x0))


    Z_mod_L = Z_modulo[:, 0:middle]
    Z_mod_R = Z_modulo[:, middle + 1:]
    X_L = X[0:middle]
    X_R = X[middle + 1:]

    N_L = integrate.simpson(Z_mod_L ** 2, X_L) / L
    N_R = integrate.simpson(Z_mod_R ** 2, X_R) / L

    N_L = filtro_array(3, N_L)
    N_R = filtro_array(3, N_R)

    np.savetxt(directory + '/T_ts.txt', T, delimiter=',')
    np.savetxt(directory + '/NL_ts.txt', N_L, delimiter=',')
    np.savetxt(directory + '/NR_ts.txt', N_R, delimiter=',')

    # legend
    fig, (ax1) = plt.subplots(1)
    pcm = ax1.pcolormesh(X - x0, T, Z_modulo, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    ax1.set_xlabel('$x$', size='25')
    ax1.tick_params(labelsize=18)
    ax1.set_ylabel('$t/T$', size='25')
    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)
    plt.tight_layout()
    plt.savefig("ST_lines.png", dpi=200)
    plt.xlim([-150, 150])
    plt.show()
    plt.close()


    fig, ax1 = plt.subplots()
    ax1.plot(T, N_R, label="$N_R(t)$", c="b", linewidth=2)
    ax1.plot(T, N_L, label="$N_L(t)$", c="r", linewidth=2)
    ax1.set_xlabel("$t/T$", fontsize=25)
    ax1.set_xlim(0, 400)
    #ax1.set_ylim(0, 6)
    ax1.set_aspect(0.4 * (400)/3)
    ax1.legend(loc="upper right", fontsize=15)
    ax1.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(directory + "/timeseries_N.png", dpi=200)
    plt.show()
    plt.close()