import matplotlib.pyplot as plt

from back_process import *


if __name__ == "__main__":
    principal_sim = "E:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras"
    principal_exp = "E:/mnustes_science/experimental_data"
    dir_01 = principal_sim + "/dimensional_dist_02/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=31.000"#/dist_03/alpha=13.0480/beta=1.000/mu=0.0950/nu=0.0180/sigma=12.6700/gamma=0.1690/dist=40.000"#
    dir_02 = principal_sim + "/dimensional_dist_02/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=21.0000"#"/dist_03/alpha=13.0480/beta=1.000/mu=0.0950/nu=0.0180/sigma=12.6700/gamma=0.1690/dist=31.000"#
    dir_03 = principal_exp + "/PT_04/f=14.20_a=10.00"
    dir_04 = principal_exp + "/PT_06/f=14.20_a=10.50_d2"
    beta = 0.004811649356064012

    ###

    Z_r_01 = np.loadtxt(dir_01 + '/field_real.txt', delimiter=',')
    Z_i_01 = np.loadtxt(dir_01 + '/field_img.txt', delimiter=',')
    X_01 = np.loadtxt(dir_01 + '/X.txt', delimiter=',')
    T_01 = np.loadtxt(dir_01 + '/T.txt', delimiter=',')
    Z_complex_01 = Z_r_01 + 1j * Z_i_01
    arg_01 = np.angle(Z_complex_01)
    Z_modulo_01 = np.absolute(Z_complex_01) / np.sqrt(beta)
    #Z_modulo_01 = filtro_superficie(Z_modulo_01, 5, "XY")
    x01 = 0

    middle = np.argmin(np.abs(X_01 - x01))
    Z_mod_L = Z_modulo_01[:, 0:middle]
    Z_mod_R = Z_modulo_01[:, middle + 1:]
    X_01_L = X_01[0:middle]
    X_01_R = X_01[middle + 1:]
    N_01_L = integrate.simpson(Z_mod_L ** 2, X_01_L) / 6
    N_01_R = integrate.simpson(Z_mod_R ** 2, X_01_R) / 6
    #N_01_L = filtro_array(3, N_01_L)
    #N_01_R = filtro_array(3, N_01_R)

    ###

    Z_r_02 = np.loadtxt(dir_02 + '/field_real.txt', delimiter=',')
    Z_i_02 = np.loadtxt(dir_02 + '/field_img.txt', delimiter=',')
    X_02 = np.loadtxt(dir_02 + '/X.txt', delimiter=',')
    T_02 = np.loadtxt(dir_02 + '/T.txt', delimiter=',')
    Z_complex_02 = Z_r_02 + 1j * Z_i_02
    arg_02 = np.angle(Z_complex_02)
    Z_modulo_02 = np.absolute(Z_complex_02) / np.sqrt(beta)
    #Z_modulo_02 = filtro_superficie(Z_modulo_02, 5, "XY")
    x02 = 0

    middle = np.argmin(np.abs(X_02 - x02))
    Z_mod_L = Z_modulo_02[:, 0:middle]
    Z_mod_R = Z_modulo_02[:, middle + 1:]
    X_02_L = X_02[0:middle]
    X_02_R = X_02[middle + 1:]
    N_02_L = integrate.simpson(Z_mod_L ** 2, X_02_L) / 6
    N_02_R = integrate.simpson(Z_mod_R ** 2, X_02_R) / 6
    #N_02_L = filtro_array(3, N_02_L)
    #N_02_R = filtro_array(3, N_02_R)

    ###

    g = 9790
    l_y = 16
    d = 20
    k_y = np.pi / l_y
    k = k_y
    tau = np.tanh(k * d)
    w_1 = np.sqrt(g * k * tau)

    Z_03 = np.loadtxt(dir_03 + '/Z_strobo.txt', delimiter=',')
    X_03 = np.loadtxt(dir_03 + '/X_mm.txt', delimiter=',')
    T_03 = np.loadtxt(dir_03 + '/T_strobo.txt', delimiter=',')
    injection_03 = np.loadtxt(dir_03 + '/IL_mm.txt', delimiter=',')
    Z_03 = filtro_superficie(Z_03, 5, "YX")
    x03 = -7
    X_03 = X_03 - x03
    T_03 = T_03 * (w_1 / (2 * np.pi))

    middle = np.argmin(np.abs(X_03 - x03))
    Z_mod_L = Z_03[:, 0:middle]
    Z_mod_R = Z_03[:, middle + 1:]
    X_03_L = X_03[0:middle]
    X_03_R = X_03[middle + 1:]
    N_03_L = integrate.simpson(Z_mod_L ** 2, X_03_L) / 13
    N_03_R = integrate.simpson(Z_mod_R ** 2, X_03_R) / 13
    #N_03_L = filtro_array(3, N_03_L)
    #N_03_R = filtro_array(3, N_03_R)

    ###

    Z_04 = np.loadtxt(dir_04 + '/Z_strobo.txt', delimiter=',')
    X_04 = np.loadtxt(dir_04 + '/X_mm.txt', delimiter=',')
    T_04 = np.loadtxt(dir_04 + '/T_strobo.txt', delimiter=',')
    injection_04 = np.loadtxt(dir_04 + '/IL_mm.txt', delimiter=',')
    Z_04 = filtro_superficie(Z_04, 5, "YX")
    x04 = injection_04[0] + injection_04[1]
    X_04 = X_04 - x04
    T_04 = T_04 * (w_1 / (2 * np.pi))

    middle = np.argmin(np.abs(X_04 - x04))
    Z_mod_L = Z_04[:, 0:middle]
    Z_mod_R = Z_04[:, middle + 1:]
    X_04_L = X_04[0:middle]
    X_04_R = X_04[middle + 1:]
    N_04_L = integrate.simpson(Z_mod_L ** 2, X_04_L) / 13
    N_04_R = integrate.simpson(Z_mod_R ** 2, X_04_R) / 13
    #N_04_L = filtro_array(3, N_04_L)
    #N_04_R = filtro_array(3, N_04_R)

    ###

    import matplotlib.gridspec as gridspec

    #fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)
    fig = plt.figure()
    outer = gridspec.GridSpec(2, 1, wspace=0.2, hspace=0.2)

    inner_01 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0], wspace=0.1, hspace=0.)
    inner_02 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[1], wspace=0.1, hspace=0.)

    ax1 = plt.Subplot(fig, inner_01[0])
    ax2 = plt.Subplot(fig, inner_01[1])
    ax3 = plt.Subplot(fig, inner_02[0])
    ax4 = plt.Subplot(fig, inner_02[1])

    norm = 500
    ax1.plot(T_03, N_03_R / (norm / 13), label="$N_R$", c="b", linewidth=2)
    ax1.plot(T_03, N_03_L / (norm / 13), label="$N_L$", c="r", linewidth=2)
    ax1.set_xlim(0, 400)
    ax1.grid(alpha=0.3)
    ax1.set_ylim(0, 1.4)
    ax1.set_xticklabels([])
    # ax1.set_yticklabels([])
    ax1.tick_params(axis="x", direction="in", top=True, bottom=False)
    ax1.tick_params(axis="y", direction="in", labelsize=20)
    xi_sim = 0
    xf_sim = 180
    ti_sim = 0
    tf_sim = 1.4
    ax1.fill_between([xi_sim, xi_sim + 0.305 * (xf_sim - xi_sim)], [0.205 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax1.plot([xi_sim, xi_sim + 0.3 * (xf_sim - xi_sim)], [0.2 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax1.plot([xi_sim, xi_sim + 0.3 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax1.plot([xi_sim, xi_sim], [0.2 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax1.text(xi_sim + 0.01 * (xf_sim - xi_sim), ti_sim + 0.43 * (tf_sim - ti_sim), "$\\textbf{EXP}$", fontsize=17,rotation=45, zorder=13)

    ax2.plot(T_04, N_04_R / (norm / 13), label="$N_R$", c="b", linewidth=2)
    ax2.plot(T_04, N_04_L / (norm / 13), label="$N_L$", c="r", linewidth=2)
    ax2.grid(alpha=0.3)
    ax2.set_xlim(0, 400)
    ax2.set_ylim(0, 1.4)
    ax2.set_xticklabels([])
    # ax2.set_yticklabels([])
    ax2.tick_params(axis="x", direction="in")
    ax2.tick_params(axis="y", direction="in", labelsize=20)

    ax3.plot(T_02 - 2050, N_02_R / (norm / 6), label="$N_R$", c="b", linewidth=2)
    ax3.plot(T_02 - 2050, N_02_L / (norm / 6), label="$N_L$", c="r", linewidth=2)
    ax3.set_xlim(0, 400)
    ax3.grid(alpha=0.3)
    ax3.set_ylim(0, 1.4)
    ax3.set_xticklabels([])
    #ax3.set_yticklabels([])
    ax3.tick_params(axis="x", direction="in", top=True, bottom=False)
    ax3.tick_params(axis="y", direction="in", labelsize=20)
    xi_sim = 0
    xf_sim = 180
    ti_sim = 0
    tf_sim = 1.4
    ax3.fill_between([xi_sim, xi_sim + 0.305 * (xf_sim - xi_sim)], [0.205 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax3.plot([xi_sim, xi_sim + 0.3 * (xf_sim - xi_sim)], [0.2 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax3.plot([xi_sim, xi_sim + 0.3 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax3.plot([xi_sim, xi_sim], [0.2 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax3.text(xi_sim + 0.008 * (xf_sim - xi_sim), ti_sim + 0.42 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=17,rotation=45, zorder=13)

    ax4.plot(T_01 - 2000, N_01_R / (norm / 6), label="$N_R$", c="b", linewidth=2)
    ax4.plot(T_01 - 2000, N_01_L / (norm / 6), label="$N_L$",  c="r", linewidth=2)
    ax4.set_xlabel("$t/T$", fontsize=20)
    ax4.set_xlim(0, 400)
    ax4.set_ylim(0, 1.4)
    ax4.grid(alpha=0.3)
    #ax4.set_yticklabels([])
    ax4.tick_params(axis="x", direction="in", labelsize=20)
    ax4.tick_params(axis="y", direction="in", labelsize=20)

    props = dict(boxstyle='round,pad=0.1', facecolor='white', alpha=0.75)
    ax1.text(362, 0.25, "$\\textbf{(a)}$", fontsize=23, rotation=0, zorder=13, bbox=props)
    ax2.text(362, 0.25, "$\\textbf{(b)}$", fontsize=23, rotation=0, zorder=13, bbox=props)
    ax3.text(362, 0.25, "$\\textbf{(c)}$", fontsize=23, rotation=0, zorder=13, bbox=props)
    ax4.text(362, 0.25, "$\\textbf{(d)}$", fontsize=23, rotation=0, zorder=13, bbox=props)


    #ax2.legend(loc="upper right", fontsize=18)
    ax4.legend(loc="center left", fontsize=16)
    fig.add_subplot(ax1)
    fig.add_subplot(ax2)
    fig.add_subplot(ax3)
    fig.add_subplot(ax4)
    fig.set_tight_layout(True)
    #plt.show()
    plt.savefig('Fig02a_alt.png', dpi=300)
    plt.close()


