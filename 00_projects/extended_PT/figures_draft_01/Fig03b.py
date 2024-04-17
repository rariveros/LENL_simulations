import matplotlib.pyplot as plt

from back_process import *


if __name__ == "__main__":
    principal_sim = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras"
    dir_01 = principal_sim + "/dimensional_dist_02/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=21.0000"#"/dist_03/alpha=13.0480/beta=1.000/mu=0.0950/nu=0.0180/sigma=12.6700/gamma=0.1690/dist=31.000"#
    beta = 0.004811649356064012

    ###

    Z_r_01 = np.loadtxt(dir_01 + '/field_real.txt', delimiter=',')
    Z_i_01 = np.loadtxt(dir_01 + '/field_img.txt', delimiter=',')
    X_01 = np.loadtxt(dir_01 + '/X.txt', delimiter=',')
    T_01 = np.loadtxt(dir_01 + '/T.txt', delimiter=',')
    Z_complex_01 = Z_r_01 + 1j * Z_i_01
    arg_01 = np.angle(Z_complex_01)
    Z_modulo_01 = np.absolute(Z_complex_01)
    Z_modulo_01 = filtro_superficie(Z_modulo_01, 5, "XY") / np.sqrt(beta)
    x01 = 0

    I_left = np.argmin(np.abs(X_01 + 21 / 2))
    I_right = np.argmin(np.abs(X_01 - 21 / 2))

    ###

    fig, ((ax1), (ax2), (ax3)) = plt.subplots(nrows=3, ncols=1)
    ax1.plot(T_01 - 2000, Z_r_01[:, I_right], label="$\\textrm{Re}\{A(x_R)\}$", c="b", linewidth=2)
    ax1.plot(T_01 - 2000, Z_i_01[:, I_right], label="$\\textrm{Im}\{A(x_R)\}$", c="b", linewidth=2, linestyle="--")
    ax1.set_xlabel("$t/T$", fontsize=25)
    ax1.set_xlim(0, 400)
    ax1.set_ylim(- 1.1 * np.amax(Z_r_01[:, I_right]), 1.1 * np.amax(Z_r_01[:, I_right]))
    ax1.grid(alpha=0.3)
    ax1.set_yticklabels([])
    ax1.tick_params(axis="x", direction="in", labelsize=15)
    ax1.tick_params(axis="y", direction="in")

    ax2.plot(T_01 - 2000, Z_r_01[:, I_left], label="$\\textrm{Re}\{A(x_L)\}$", c="r", linewidth=2)
    ax2.plot(T_01 - 2000, Z_i_01[:, I_left], label="$\\textrm{Im}\{A(x_L)\}$", c="r", linewidth=2, linestyle="--")
    ax2.set_xlabel("$t/T$", fontsize=25)
    ax2.set_xlim(0, 400)
    ax2.set_ylim(- 1.1 * np.amax(Z_r_01[:, I_right]), 1.1 * np.amax(Z_r_01[:, I_right]))
    ax2.grid(alpha=0.3)
    ax2.set_yticklabels([])
    ax2.tick_params(axis="x", direction="in", labelsize=15)
    ax2.tick_params(axis="y", direction="in")

    ax1.legend(loc="upper right", fontsize=13)
    ax2.legend(loc="upper right", fontsize=13)
    plt.tight_layout()
    #plt.show()
    plt.savefig('Fig03b.png', dpi=300)
    plt.close()