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
    x0 = 0
    inj_L = injection[1] - x0
    inj_R = injection[0] - x0

    d = "dist_02"

    if d == "dist_01":
        medio = (inj_L + inj_R) - 22 / 2
        dist_R = X - (medio + 16)
        dist_L = X - (medio - 16)
    elif d == "dist_02":
        medio = ((inj_L + inj_R) / 2)
        dist_R = X - (medio + 16)
        dist_L = X - (medio - 16)
    J_L = np.argmin(np.abs(dist_L))
    J_R = np.argmin(np.abs(dist_R))

    Z_modulo = filtro_superficie(Z_modulo, 2, "Y")

    np.savetxt(directory+ '/T_ts.txt', T, delimiter=',')
    np.savetxt(directory + '/AL_ts.txt', Z_modulo[:, J_L], delimiter=',')
    np.savetxt(directory + '/AR_ts.txt', Z_modulo[:, J_R], delimiter=',')

    matplotlib.rc('xtick', labelsize=15)
    matplotlib.rc('ytick', labelsize=15)
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
    ax1.vlines([X[J_L], X[J_R]], 0, T[-1], colors=["r", "b"], linestyles="--")
    plt.savefig("ST_lines.png", dpi=200)
    plt.xlim([-150, 150])
    ax1.vlines([inj_L, inj_R, medio], 0, T[-1], colors="g", linestyles="-")
    plt.show()
    plt.close()


    fig, (ax1, ax2, ax3) = plt.subplots(3)
    ax1.plot(T, Z_modulo[:, J_R], label="$|A(x_R, t)|$", c="b", linewidth=2)
    ax1.plot(T, Z_modulo[:, J_L], label="$|A(x_L, t)|$", c="r", linewidth=2)
    ax1.set_xlabel("$t/T$", fontsize=25)
    ax1.set_xlim(0, 400)
    ax1.set_ylim(0, 5.1)
    #ax1.set_aspect(0.4 * (400)/4.2)
    ax1.legend(loc="upper right", fontsize=15)
    ax1.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(directory + "/timeseries_A.png", dpi=300)
    plt.show()
    plt.close()
