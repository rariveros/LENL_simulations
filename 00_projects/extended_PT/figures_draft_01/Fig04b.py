import matplotlib.pyplot as plt

from back_process import *


if __name__ == "__main__":
    principal_sim = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras"
    principal_exp = "C:/mnustes_science/experimental_data"
    dir_05 = principal_sim + "/dimensional/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.2800/dist=20.000"
    dir_06 = principal_sim + "/dimensional/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.2200/dist=20.000"
    dir_07 = principal_exp + "/PT_05/f=14.20_a=12.00"
    dir_08 = principal_exp + "/PT_05/f=14.20_a=14.00"

    beta = 0.004811649356064012

    Z_r_05 = np.loadtxt(dir_05 + '/field_real.txt', delimiter=',')
    Z_i_05 = np.loadtxt(dir_05 + '/field_img.txt', delimiter=',')
    X_05 = np.loadtxt(dir_05 + '/X.txt', delimiter=',')
    T_05 = np.loadtxt(dir_05 + '/T.txt', delimiter=',')
    Z_complex_05 = Z_r_05 + 1j * Z_i_05
    arg_05 = np.angle(Z_complex_05)
    Z_modulo_05 = np.absolute(Z_complex_05)
    Z_modulo_05 = filtro_superficie(Z_modulo_05, 2, "XY") / np.sqrt(beta)

    x05 = 0
    middle = np.argmin(np.abs(X_05 - x05))
    Z_mod_L = Z_modulo_05[:, 0:middle]
    Z_mod_R = Z_modulo_05[:, middle + 1:]
    X_05_L = X_05[0:middle]
    X_05_R = X_05[middle + 1:]
    N_05_L = integrate.simpson(Z_mod_L ** 2, X_05_L)
    N_05_R = integrate.simpson(Z_mod_R ** 2, X_05_R)
    N_05_L = filtro_array(3, N_05_L)
    N_05_R = filtro_array(3, N_05_R)

    ###

    Z_r_06 = np.loadtxt(dir_06 + '/field_real.txt', delimiter=',')
    Z_i_06 = np.loadtxt(dir_06 + '/field_img.txt', delimiter=',')
    X_06 = np.loadtxt(dir_06 + '/X.txt', delimiter=',')
    T_06 = np.loadtxt(dir_06 + '/T.txt', delimiter=',')
    Z_complex_06 = Z_r_06 + 1j * Z_i_06
    arg_06 = np.angle(Z_complex_06)
    Z_modulo_06 = np.absolute(Z_complex_06)
    Z_modulo_06 = filtro_superficie(Z_modulo_06, 2, "XY") / np.sqrt(beta)
    x06 = 0

    middle = np.argmin(np.abs(X_06 - x06))
    Z_mod_L = Z_modulo_06[:, 0:middle]
    Z_mod_R = Z_modulo_06[:, middle + 1:]
    X_06_L = X_06[0:middle]
    X_06_R = X_06[middle + 1:]
    N_06_L = integrate.simpson(Z_mod_L ** 2, X_06_L)
    N_06_R = integrate.simpson(Z_mod_R ** 2, X_06_R)
    N_06_L = filtro_array(3, N_06_L)
    N_06_R = filtro_array(3, N_06_R)

    ###

    g = 9790
    l_y = 16
    d = 20
    k_y = np.pi / l_y
    k = k_y
    tau = np.tanh(k * d)
    w_1 = np.sqrt(g * k * tau)

    Z_07 = np.loadtxt(dir_07 + '/Z_strobo.txt', delimiter=',')
    X_07 = np.loadtxt(dir_07 + '/X_mm.txt', delimiter=',')
    T_07 = np.loadtxt(dir_07 + '/T_strobo.txt', delimiter=',')
    injection_07 = np.loadtxt(dir_07 + '/IL_mm.txt', delimiter=',')
    #Z_07 = filtro_superficie(Z_07, 3, "YX")

    x07 = injection_07[0] + injection_07[1] - 16 * 2.5
    X_07 = X_07 - x07
    T_07 = T_07 * (w_1 / (2 * np.pi))

    middle = np.argmin(np.abs(X_07 - x07))
    Z_mod_L = Z_07[:, 0:middle]
    Z_mod_R = Z_07[:, middle + 1:]
    X_07_L = X_07[0:middle]
    X_07_R = X_07[middle + 1:]
    N_07_L = integrate.simpson(Z_mod_L ** 2, X_07_L)
    N_07_R = integrate.simpson(Z_mod_R ** 2, X_07_R)
    N_07_L = filtro_array(3, N_07_L)
    N_07_R = filtro_array(3, N_07_R)

    ###

    Z_08 = np.loadtxt(dir_08 + '/Z_strobo.txt', delimiter=',')
    X_08 = np.loadtxt(dir_08 + '/X_mm.txt', delimiter=',')
    T_08 = np.loadtxt(dir_08 + '/T_strobo.txt', delimiter=',')
    injection_08 = np.loadtxt(dir_08 + '/IL_mm.txt', delimiter=',')
    #Z_08 = filtro_superficie(Z_08, 3, "YX")

    x08 = injection_08[0] + injection_08[1] - 16 * 3 + 0
    X_08 = X_08 - x08
    T_08 = T_08 * (w_1 / (2 * np.pi)) * 2

    middle = np.argmin(np.abs(X_08 - x08))
    Z_mod_L = Z_08[:, 0:middle]
    Z_mod_R = Z_08[:, middle + 1:]
    X_08_L = X_08[0:middle]
    X_08_R = X_08[middle + 1:]
    N_08_L = integrate.simpson(Z_mod_L ** 2, X_08_L)
    N_08_R = integrate.simpson(Z_mod_R ** 2, X_08_R)
    N_08_L = filtro_array(3, N_08_L)
    N_08_R = filtro_array(3, N_08_R)

    ###

    fig, (ax2, ax1) = plt.subplots(2)

    ax1.plot(T_06 - 2103, N_06_R, label="$N_R^{(\\textrm{N})}$", c="b", linestyle="--", linewidth=2)
    ax1.plot(T_06 - 2103, N_06_L, label="$N_L^{(\\textrm{N})}$", c="r", linestyle="--", linewidth=2)
    ax1.plot(T_07, N_07_R, label="$N_R^{(\\textrm{E})}$", c="b", linewidth=2)
    ax1.plot(T_07, N_07_L, label="$N_L^{(\\textrm{E})}$", c="r", linewidth=2)
    ax1.set_xlim(0, 200)
    ax1.grid(alpha=0.3)
    ax1.set_xlabel("$t/T$", fontsize=25)
    ax1.set_ylabel("$N/N_{\\textrm{max}}$", fontsize=25)
    ax1.set_yticklabels([])
    ax1.tick_params(axis="x", direction="in", labelsize=15)
    ax1.tick_params(axis="y", direction="in")
    xi_sim = 0
    xf_sim = 180
    ti_sim = 0
    tf_sim = 1.1 * np.amax(N_06_R)

    ax2.plot(T_05 - 2015, N_05_R, label="$N_R\ (\\textrm{NUM})$", c="b", linestyle="--", linewidth=2)
    ax2.plot(T_05 - 2015, N_05_L, label="$N_L\ (\\textrm{NUM})$", c="r", linestyle="--", linewidth=2)
    ax2.plot(T_08, N_08_R, label="$N_R(t)$", c="b", linewidth=2)
    ax2.plot(T_08, N_08_L, label="$N_L(t)$", c="r", linewidth=2)
    ax2.grid(alpha=0.3)
    ax2.set_xlim(0, 200)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.tick_params(axis="x", direction="in")
    ax2.tick_params(axis="y", direction="in")

    ax1.legend(loc="upper right", fontsize=13)
    plt.tight_layout()
    #plt.show()
    plt.savefig('Fig04b.png', dpi=300)
    plt.close()


